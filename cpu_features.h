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

/*
 * We safely assume every device that we support supports crc32/pmull hw
 * features. Practically speaking, we never supported devices without it.
 *
 * Hence, simply enable ARMV8_ALWAYS_NEW to avoid the cost of evaluating
 * it during runtime.
 *
 * In the 0% chance of us adding support for devices without these hw
 * features, it'll be very apparent during bring-up and not cause a
 * silent breakage.
 *
 * If you encounter build breakage with this change, you should review
 * your compiler flag and see if it correctly enables crc hw feature.
 */
#define ARMV8_ALWAYS_NEW

#if defined(ARMV8_ALWAYS_NEW)
#define arm_cpu_enable_crc32 1
#define arm_cpu_enable_pmull 1
#define cpu_check_features() ((void)0)
#else
extern int arm_cpu_enable_crc32;
extern int arm_cpu_enable_pmull;
void cpu_check_features(void);
#endif
