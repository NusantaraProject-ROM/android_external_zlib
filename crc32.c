/* crc32.c -- compute the CRC-32 of a data stream
 * Copyright (C) 1995-2006, 2010, 2011, 2012, 2016, 2018 Mark Adler
 * For conditions of distribution and use, see copyright notice in zlib.h
 *
 * Thanks to Rodney Brown <rbrown64@csc.com.au> for his contribution of faster
 * CRC methods: exclusive-oring 32 bits of data at a time, and pre-computing
 * tables for updating the shift register in one step with three exclusive-ors
 * instead of four steps with four exclusive-ors.  This results in about a
 * factor of two increase in speed on a Power PC G4 (PPC7455) using gcc -O3.
 */

/* @(#) $Id$ */

/*
  Note on the use of DYNAMIC_CRC_TABLE: there is no mutex or semaphore
  protection on the static variables used to control the first-use generation
  of the crc tables.  Therefore, if you #define DYNAMIC_CRC_TABLE, you should
  first call get_crc_table() to initialize the tables before allowing more than
  one thread to use crc32().

  DYNAMIC_CRC_TABLE and MAKECRCH can be #defined to write out crc32.h. A main()
  routine is also produced, so that this one source file can be compiled to an
  executable.
 */

#ifdef MAKECRCH
#  include <stdio.h>
#  ifndef DYNAMIC_CRC_TABLE
#    define DYNAMIC_CRC_TABLE
#  endif /* !DYNAMIC_CRC_TABLE */
#endif /* MAKECRCH */

#include "deflate.h"
#include "x86.h"
#include "zutil.h"      /* for STDC and FAR definitions */

#if defined(CRC32_SIMD_SSE42_PCLMUL)
#include "crc32_simd.h"
#elif defined(CRC32_ARMV8_CRC32)
#include "arm_features.h"
#include "crc32_simd.h"
#endif

/* Definitions for doing the crc four data bytes at a time. */
#if !defined(NOBYFOUR) && defined(Z_U4)
#  define BYFOUR
#endif
#ifdef BYFOUR
   local unsigned long crc32_little OF((unsigned long,
                        const unsigned char FAR *, z_size_t));
   local unsigned long crc32_big OF((unsigned long,
                        const unsigned char FAR *, z_size_t));
#  define TBLS 8
#else
#  define TBLS 1
#endif /* BYFOUR */

/* Local functions for crc concatenation */
#define GF2_DIM 32      /* dimension of GF(2) vectors (length of CRC) */
local z_crc_t gf2_matrix_times OF((const z_crc_t *mat, z_crc_t vec));
local uLong crc32_combine_ OF((uLong crc1, uLong crc2, z_off64_t len2));

/* ========================================================================= */
local z_crc_t gf2_matrix_times(mat, vec)
    const z_crc_t *mat;
    z_crc_t vec;
{
    z_crc_t sum;

    sum = 0;
    while (vec) {
        if (vec & 1)
            sum ^= *mat;
        vec >>= 1;
        mat++;
    }
    return sum;
}


#ifdef DYNAMIC_CRC_TABLE

local volatile int crc_table_empty = 1;
local z_crc_t FAR crc_table[TBLS][256];
local z_crc_t FAR crc_comb[GF2_DIM][GF2_DIM];
local void make_crc_table OF((void));
local void gf2_matrix_square OF((z_crc_t *square, const z_crc_t *mat));
#ifdef MAKECRCH
   local void write_table OF((FILE *, const z_crc_t FAR *, int));
#endif /* MAKECRCH */

/* ========================================================================= */
local void gf2_matrix_square(square, mat)
    z_crc_t *square;
    const z_crc_t *mat;
{
    int n;

    for (n = 0; n < GF2_DIM; n++)
        square[n] = gf2_matrix_times(mat, mat[n]);
}

/*
  Generate tables for a byte-wise 32-bit CRC calculation on the polynomial:
  x^32+x^26+x^23+x^22+x^16+x^12+x^11+x^10+x^8+x^7+x^5+x^4+x^2+x+1.

  Polynomials over GF(2) are represented in binary, one bit per coefficient,
  with the lowest powers in the most significant bit.  Then adding polynomials
  is just exclusive-or, and multiplying a polynomial by x is a right shift by
  one.  If we call the above polynomial p, and represent a byte as the
  polynomial q, also with the lowest power in the most significant bit (so the
  byte 0xb1 is the polynomial x^7+x^3+x+1), then the CRC is (q*x^32) mod p,
  where a mod b means the remainder after dividing a by b.

  This calculation is done using the shift-register method of multiplying and
  taking the remainder.  The register is initialized to zero, and for each
  incoming bit, x^32 is added mod p to the register if the bit is a one (where
  x^32 mod p is p+x^32 = x^26+...+1), and the register is multiplied mod p by
  x (which is shifting right by one and adding x^32 mod p if the bit shifted
  out is a one).  We start with the highest power (least significant bit) of
  q and repeat for all eight bits of q.

  The first table is simply the CRC of all possible eight bit values.  This is
  all the information needed to generate CRCs on data a byte at a time for all
  combinations of CRC register values and incoming bytes.  The remaining tables
  allow for word-at-a-time CRC calculation for both big-endian and little-
  endian machines, where a word is four bytes.
*/
local void make_crc_table()
{
    z_crc_t c;
    int n, k;
    z_crc_t poly;                       /* polynomial exclusive-or pattern */
    /* terms of polynomial defining this crc (except x^32): */
    static volatile int first = 1;      /* flag to limit concurrent making */
    static const unsigned char p[] = {0,1,2,4,5,7,8,10,11,12,16,22,23,26};

    /* See if another task is already doing this (not thread-safe, but better
       than nothing -- significantly reduces duration of vulnerability in
       case the advice about DYNAMIC_CRC_TABLE is ignored) */
    if (first) {
        first = 0;

        /* make exclusive-or pattern from polynomial (0xedb88320UL) */
        poly = 0;
        for (n = 0; n < (int)(sizeof(p)/sizeof(unsigned char)); n++)
            poly |= (z_crc_t)1 << (31 - p[n]);

        /* generate a crc for every 8-bit value */
        for (n = 0; n < 256; n++) {
            c = (z_crc_t)n;
            for (k = 0; k < 8; k++)
                c = c & 1 ? poly ^ (c >> 1) : c >> 1;
            crc_table[0][n] = c;
        }

#ifdef BYFOUR
        /* generate crc for each value followed by one, two, and three zeros,
           and then the byte reversal of those as well as the first table */
        for (n = 0; n < 256; n++) {
            c = crc_table[0][n];
            crc_table[4][n] = ZSWAP32(c);
            for (k = 1; k < 4; k++) {
                c = crc_table[0][c & 0xff] ^ (c >> 8);
                crc_table[k][n] = c;
                crc_table[k + 4][n] = ZSWAP32(c);
            }
        }
#endif /* BYFOUR */

        /* generate zero operators table for crc32_combine() */

        /* generate the operator to apply a single zero bit to a CRC -- the
           first row adds the polynomial if the low bit is a 1, and the
           remaining rows shift the CRC right one bit */
        k = GF2_DIM - 3;
        crc_comb[k][0] = 0xedb88320UL;      /* CRC-32 polynomial */
        z_crc_t row = 1;
        for (n = 1; n < GF2_DIM; n++) {
            crc_comb[k][n] = row;
            row <<= 1;
        }

        /* generate operators that apply 2, 4, and 8 zeros to a CRC, putting
           the last one, the operator for one zero byte, at the 0 position */
        gf2_matrix_square(crc_comb[k + 1], crc_comb[k]);
        gf2_matrix_square(crc_comb[k + 2], crc_comb[k + 1]);
        gf2_matrix_square(crc_comb[0], crc_comb[k + 2]);

        /* generate operators for applying 2^n zero bytes to a CRC, filling out
           the remainder of the table -- the operators repeat after GF2_DIM
           values of n, so the table only needs GF2_DIM entries, regardless of
           the size of the length being processed */
        for (n = 1; n < k; n++)
            gf2_matrix_square(crc_comb[n], crc_comb[n - 1]);

        /* mark tables as complete, in case someone else is waiting */
        crc_table_empty = 0;
    }
    else {      /* not first */
        /* wait for the other guy to finish (not efficient, but rare) */
        while (crc_table_empty)
            ;
    }
#ifdef MAKECRCH
    {
        FILE *out;

        out = fopen("crc32.h", "w");
        if (out == NULL) return;

        /* write out CRC table to crc32.h */
        fprintf(out, "/* crc32.h -- tables for rapid CRC calculation\n");
        fprintf(out, " * Generated automatically by crc32.c\n */\n\n");
        fprintf(out, "local const z_crc_t FAR ");
        fprintf(out, "crc_table[%d][256] =\n{\n  {\n", TBLS);
        write_table(out, crc_table[0], 256);
#  ifdef BYFOUR
        fprintf(out, "#ifdef BYFOUR\n");
        for (k = 1; k < 8; k++) {
            fprintf(out, "  },\n  {\n");
            write_table(out, crc_table[k], 256);
        }
        fprintf(out, "#endif\n");
#  endif /* BYFOUR */
        fprintf(out, "  }\n};\n");

        /* write out zero operator table to crc32.h */
        fprintf(out, "\nlocal const z_crc_t FAR ");
        fprintf(out, "crc_comb[%d][%d] =\n{\n  {\n", GF2_DIM, GF2_DIM);
        write_table(out, crc_comb[0], GF2_DIM);
        for (k = 1; k < GF2_DIM; k++) {
            fprintf(out, "  },\n  {\n");
            write_table(out, crc_comb[k], GF2_DIM);
        }
        fprintf(out, "  }\n};\n");
        fclose(out);
    }
#endif /* MAKECRCH */
}

#ifdef MAKECRCH
local void write_table(out, table, k)
    FILE *out;
    const z_crc_t FAR *table;
    int k;
{
    int n;

    for (n = 0; n < k; n++)
        fprintf(out, "%s0x%08lxUL%s", n % 5 ? "" : "    ",
                (unsigned long)(table[n]),
                n == k - 1 ? "\n" : (n % 5 == 4 ? ",\n" : ", "));
}

int main()
{
    make_crc_table();
    return 0;
}
#endif /* MAKECRCH */

#else /* !DYNAMIC_CRC_TABLE */
/* ========================================================================
 * Tables of CRC-32s of all single-byte values, made by make_crc_table(),
 * and tables of zero operator matrices for crc32_combine().
 */
#include "crc32.h"
#endif /* DYNAMIC_CRC_TABLE */

/* =========================================================================
 * This function can be used by asm versions of crc32()
 */
const z_crc_t FAR * ZEXPORT get_crc_table()
{
#ifdef DYNAMIC_CRC_TABLE
    if (crc_table_empty)
        make_crc_table();
#endif /* DYNAMIC_CRC_TABLE */
    return (const z_crc_t FAR *)crc_table;
}

/* ========================================================================= */
#define DO1 crc = crc_table[0][((int)crc ^ (*buf++)) & 0xff] ^ (crc >> 8)
#define DO8 DO1; DO1; DO1; DO1; DO1; DO1; DO1; DO1

/* ========================================================================= */

/* =========================================================================
 * Use ARM machine instructions if available. This will compute the CRC about
 * ten times faster than the braided calculation. This code does not check for
 * the presence of the CRC instruction at run time. __ARM_FEATURE_CRC32 will
 * only be defined if the compilation specifies an ARM processor architecture
 * that has the instructions. For example, compiling with -march=armv8.1-a or
 * -march=armv8-a+crc, or -march=native if the compile machine has the crc32
 * instructions.
 */
#if defined(__aarch64__) && defined(__ARM_FEATURE_CRC32) && W == 8

/*
   Constants empirically determined to maximize speed. These values are from
   measurements on a Cortex-A57. Your mileage may vary.
 */
#define Z_BATCH 3990                /* number of words in a batch */
#define Z_BATCH_ZEROS 0xa10d3d0c    /* computed from Z_BATCH = 3990 */
#define Z_BATCH_MIN 800             /* fewest words in a final batch */

unsigned long ZEXPORT crc32_z(crc, buf, len)
    unsigned long crc;
    const unsigned char FAR *buf;
    z_size_t len;
{
    z_crc_t val;
    z_word_t crc1, crc2;
    const z_word_t *word;
    z_word_t val0, val1, val2;
    z_size_t last, last2, i;
    z_size_t num;

    /* Return initial CRC, if requested. */
    if (buf == Z_NULL) return 0;

#ifdef DYNAMIC_CRC_TABLE
    once(&made, make_crc_table);
#endif /* DYNAMIC_CRC_TABLE */

    /* Pre-condition the CRC */
    crc ^= 0xffffffff;

    /* Compute the CRC up to a word boundary. */
    while (len && ((z_size_t)buf & 7) != 0) {
        len--;
        val = *buf++;
        __asm__ volatile("crc32b %w0, %w0, %w1" : "+r"(crc) : "r"(val));
    }

    /* Prepare to compute the CRC on full 64-bit words word[0..num-1]. */
    word = (z_word_t const *)buf;
    num = len >> 3;
    len &= 7;

    /* Do three interleaved CRCs to realize the throughput of one crc32x
       instruction per cycle. Each CRC is calcuated on Z_BATCH words. The three
       CRCs are combined into a single CRC after each set of batches. */
    while (num >= 3 * Z_BATCH) {
        crc1 = 0;
        crc2 = 0;
        for (i = 0; i < Z_BATCH; i++) {
            val0 = word[i];
            val1 = word[i + Z_BATCH];
            val2 = word[i + 2 * Z_BATCH];
            __asm__ volatile("crc32x %w0, %w0, %x1" : "+r"(crc) : "r"(val0));
            __asm__ volatile("crc32x %w0, %w0, %x1" : "+r"(crc1) : "r"(val1));
            __asm__ volatile("crc32x %w0, %w0, %x1" : "+r"(crc2) : "r"(val2));
        }
        word += 3 * Z_BATCH;
        num -= 3 * Z_BATCH;
        crc = multmodp(Z_BATCH_ZEROS, crc) ^ crc1;
        crc = multmodp(Z_BATCH_ZEROS, crc) ^ crc2;
    }

    /* Do one last smaller batch with the remaining words, if there are enough
       to pay for the combination of CRCs. */
    last = num / 3;
    if (last >= Z_BATCH_MIN) {
        last2 = last << 1;
        crc1 = 0;
        crc2 = 0;
        for (i = 0; i < last; i++) {
            val0 = word[i];
            val1 = word[i + last];
            val2 = word[i + last2];
            __asm__ volatile("crc32x %w0, %w0, %x1" : "+r"(crc) : "r"(val0));
            __asm__ volatile("crc32x %w0, %w0, %x1" : "+r"(crc1) : "r"(val1));
            __asm__ volatile("crc32x %w0, %w0, %x1" : "+r"(crc2) : "r"(val2));
        }
        word += 3 * last;
        num -= 3 * last;
        val = x2nmodp(last, 6);
        crc = multmodp(val, crc) ^ crc1;
        crc = multmodp(val, crc) ^ crc2;
    }

    /* Compute the CRC on any remaining words. */
    for (i = 0; i < num; i++) {
        val0 = word[i];
        __asm__ volatile("crc32x %w0, %w0, %x1" : "+r"(crc) : "r"(val0));
    }
    word += num;

    /* Complete the CRC on any remaining bytes. */
    buf = (const unsigned char FAR *)word;
    while (len) {
        len--;
        val = *buf++;
        __asm__ volatile("crc32b %w0, %w0, %w1" : "+r"(crc) : "r"(val));
    }

    /* Return the CRC, post-conditioned. */
    return crc ^ 0xffffffff;
}

#else

unsigned long ZEXPORT crc32_z(crc, buf, len)
    unsigned long crc;
    const unsigned char FAR *buf;
    z_size_t len;
{
    /*
     * zlib convention is to call crc32(0, NULL, 0); before making
     * calls to crc32(). So this is a good, early (and infrequent)
     * place to cache CPU features if needed for those later, more
     * interesting crc32() calls.
     */
#if defined(CRC32_SIMD_SSE42_PCLMUL)
    /*
     * Use x86 sse4.2+pclmul SIMD to compute the crc32. Since this
     * routine can be freely used, check CPU features here.
     */
    if (buf == Z_NULL) {
        if (!len) /* Assume user is calling crc32(0, NULL, 0); */
            x86_check_features();
        return 0UL;
    }

    if (x86_cpu_enable_simd && len >= Z_CRC32_SSE42_MINIMUM_LENGTH) {
        /* crc32 16-byte chunks */
        z_size_t chunk_size = len & ~Z_CRC32_SSE42_CHUNKSIZE_MASK;
        crc = ~crc32_sse42_simd_(buf, chunk_size, ~(uint32_t)crc);
        /* check remaining data */
        len -= chunk_size;
        if (!len)
            return crc;
        /* Fall into the default crc32 for the remaining data. */
        buf += chunk_size;
    }
#else
    if (buf == Z_NULL) {
        return 0UL;
    }
#endif /* CRC32_SIMD_SSE42_PCLMUL */

#ifdef DYNAMIC_CRC_TABLE
    if (crc_table_empty)
        make_crc_table();
#endif /* DYNAMIC_CRC_TABLE */

#ifdef BYFOUR
    if (sizeof(void *) == sizeof(ptrdiff_t)) {
        z_crc_t endian;

        endian = 1;
        if (*((unsigned char *)(&endian)))
            return crc32_little(crc, buf, len);
        else
            return crc32_big(crc, buf, len);
    }
#endif /* BYFOUR */
    crc = crc ^ 0xffffffffUL;
    while (len >= 8) {
        DO8;
        len -= 8;
    }
    if (len) do {
        DO1;
    } while (--len);
    return crc ^ 0xffffffffUL;
}
#endif

/* ========================================================================= */
unsigned long ZEXPORT crc32(crc, buf, len)
    unsigned long crc;
    const unsigned char FAR *buf;
    uInt len;
{
#if defined(CRC32_ARMV8_CRC32)
    /* We got to verify ARM CPU features, so exploit the common usage pattern
     * of calling this function with Z_NULL for an initial valid crc value.
     * This allows to cache the result of the feature check and avoid extraneous
     * function calls.
     * TODO: try to move this to crc32_z if we don't loose performance on ARM.
     */
    if (buf == Z_NULL) {
        if (!len) /* Assume user is calling crc32(0, NULL, 0); */
            arm_check_features();
        return 0UL;
    }

    if (arm_cpu_enable_crc32)
        return armv8_crc32_little(crc, buf, len);
#endif
    return crc32_z(crc, buf, len);
}

#ifdef BYFOUR

/*
   This BYFOUR code accesses the passed unsigned char * buffer with a 32-bit
   integer pointer type. This violates the strict aliasing rule, where a
   compiler can assume, for optimization purposes, that two pointers to
   fundamentally different types won't ever point to the same memory. This can
   manifest as a problem only if one of the pointers is written to. This code
   only reads from those pointers. So long as this code remains isolated in
   this compilation unit, there won't be a problem. For this reason, this code
   should not be copied and pasted into a compilation unit in which other code
   writes to the buffer that is passed to these routines.
 */

/* ========================================================================= */
#define DOLIT4 c ^= *buf4++; \
        c = crc_table[3][c & 0xff] ^ crc_table[2][(c >> 8) & 0xff] ^ \
            crc_table[1][(c >> 16) & 0xff] ^ crc_table[0][c >> 24]
#define DOLIT32 DOLIT4; DOLIT4; DOLIT4; DOLIT4; DOLIT4; DOLIT4; DOLIT4; DOLIT4

/* ========================================================================= */
local unsigned long crc32_little(crc, buf, len)
    unsigned long crc;
    const unsigned char FAR *buf;
    z_size_t len;
{
    register z_crc_t c;
    register const z_crc_t FAR *buf4;

    c = (z_crc_t)crc;
    c = ~c;
    while (len && ((ptrdiff_t)buf & 3)) {
        c = crc_table[0][(c ^ *buf++) & 0xff] ^ (c >> 8);
        len--;
    }

    buf4 = (const z_crc_t FAR *)(const void FAR *)buf;
    while (len >= 32) {
        DOLIT32;
        len -= 32;
    }
    while (len >= 4) {
        DOLIT4;
        len -= 4;
    }
    buf = (const unsigned char FAR *)buf4;

    if (len) do {
        c = crc_table[0][(c ^ *buf++) & 0xff] ^ (c >> 8);
    } while (--len);
    c = ~c;
    return (unsigned long)c;
}

/* ========================================================================= */
#define DOBIG4 c ^= *buf4++; \
        c = crc_table[4][c & 0xff] ^ crc_table[5][(c >> 8) & 0xff] ^ \
            crc_table[6][(c >> 16) & 0xff] ^ crc_table[7][c >> 24]
#define DOBIG32 DOBIG4; DOBIG4; DOBIG4; DOBIG4; DOBIG4; DOBIG4; DOBIG4; DOBIG4

/* ========================================================================= */
local unsigned long crc32_big(crc, buf, len)
    unsigned long crc;
    const unsigned char FAR *buf;
    z_size_t len;
{
    register z_crc_t c;
    register const z_crc_t FAR *buf4;

    c = ZSWAP32((z_crc_t)crc);
    c = ~c;
    while (len && ((ptrdiff_t)buf & 3)) {
        c = crc_table[4][(c >> 24) ^ *buf++] ^ (c << 8);
        len--;
    }

    buf4 = (const z_crc_t FAR *)(const void FAR *)buf;
    while (len >= 32) {
        DOBIG32;
        len -= 32;
    }
    while (len >= 4) {
        DOBIG4;
        len -= 4;
    }
    buf = (const unsigned char FAR *)buf4;

    if (len) do {
        c = crc_table[4][(c >> 24) ^ *buf++] ^ (c << 8);
    } while (--len);
    c = ~c;
    return (unsigned long)(ZSWAP32(c));
}

#endif /* BYFOUR */

/* ========================================================================= */
local uLong crc32_combine_(crc1, crc2, len2)
    uLong crc1;
    uLong crc2;
    z_off64_t len2;
{
    int n;

#ifdef DYNAMIC_CRC_TABLE
    if (crc_table_empty)
        make_crc_table();
#endif /* DYNAMIC_CRC_TABLE */

    if (len2 > 0)
        /* operator for 2^n zeros repeats every GF2_DIM n values */
        for (n = 0; len2; n = (n + 1) % GF2_DIM, len2 >>= 1)
            if (len2 & 1)
                crc1 = gf2_matrix_times(crc_comb[n], crc1);
    return crc1 ^ crc2;
}

/* ========================================================================= */
uLong ZEXPORT crc32_combine(crc1, crc2, len2)
    uLong crc1;
    uLong crc2;
    z_off_t len2;
{
    return crc32_combine_(crc1, crc2, len2);
}

uLong ZEXPORT crc32_combine64(crc1, crc2, len2)
    uLong crc1;
    uLong crc2;
    z_off64_t len2;
{
    return crc32_combine_(crc1, crc2, len2);
}

ZLIB_INTERNAL void crc_reset(deflate_state *const s)
{
    if (x86_cpu_enable_simd) {
        crc_fold_init(s);
        return;
    }
    s->strm->adler = crc32(0L, Z_NULL, 0);
}

ZLIB_INTERNAL void crc_finalize(deflate_state *const s)
{
    if (x86_cpu_enable_simd)
        s->strm->adler = crc_fold_512to32(s);
}

ZLIB_INTERNAL void copy_with_crc(z_streamp strm, Bytef *dst, long size)
{
    if (x86_cpu_enable_simd) {
        crc_fold_copy(strm->state, dst, strm->next_in, size);
        return;
    }
    zmemcpy(dst, strm->next_in, size);
    strm->adler = crc32(strm->adler, dst, size);
}
