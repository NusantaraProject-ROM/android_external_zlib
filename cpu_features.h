/* cpu_features.h -- Processor features detection.
 *
 * Copyright 2018 The Chromium Authors. All rights reserved.
 * Use of this source code is governed by a BSD-style license that can be
 * found in the Chromium source repository LICENSE file.
 */

#include "zlib.h"

/* TODO(cavalcantii): remove checks for x86_flags on deflate.
 */
extern int x86_cpu_enable_sse2;
extern int x86_cpu_enable_ssse3;
extern int x86_cpu_enable_simd;

#if defined(ARMV8_OS_MACOS)
/* Crypto extensions (crc32/pmull) are a baseline feature in ARMv8.1-A, and
 * OSX running on arm64 is new enough that these can be assumed without
 * runtime detection.
 */
#define ARMV8_ALWAYS_NEW
#endif

#if defined(ARMV8_ALWAYS_NEW)
#define arm_cpu_enable_crc32 1
#define arm_cpu_enable_pmull 1
#define cpu_check_features() ((void)0)
#else
extern int arm_cpu_enable_crc32;
extern int arm_cpu_enable_pmull;
void cpu_check_features(void);
#endif
