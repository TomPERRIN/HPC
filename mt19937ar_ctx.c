/* ------------------------------------------------------------------- */
/* mt19937ar_ctx.c                                                     */
/*                                                                     */
/* Thread-safe (context based) version of MT19937.                     */
/*                                                                     */
/* A C-program for MT19937, with initialization improved 2002/1/26.    */
/* Coded by Takuji Nishimura and Makoto Matsumoto.                     */
/*                                                                     */
/* Copyright (C) 1997 - 2002, Makoto Matsumoto and Takuji Nishimura,   */
/* All rights reserved.                                                */
/*                                                                     */
/* Redistribution and use in source and binary forms, with or without  */
/* modification, are permitted provided that the following conditions  */
/* are met:                                                            */
/*   1. Redistributions of source code must retain the above copyright */
/*      notice.                                                        */
/*   2. Redistributions in binary form must reproduce the above        */
/*      copyright notice in the documentation and/or other materials.  */
/*   3. The names of its contributors may not be used to endorse or    */
/*      promote products derived from this software without specific   */
/*      prior written permission.                                      */
/*                                                                     */
/* THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS */
/* "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT   */
/* LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS   */
/* FOR A PARTICULAR PURPOSE ARE DISCLAIMED.                            */
/*                                                                     */
/* Only difference with the original mt19937ar.c: the state vector is  */
/* a context given as a parameter instead of a global variable.        */
/* ------------------------------------------------------------------- */

#include "mt19937ar_ctx.h"

/* Period parameters */
#define N 624
#define M 397
#define MATRIX_A 0x9908b0dfUL   /* constant vector a */
#define UPPER_MASK 0x80000000UL /* most significant w-r bits */
#define LOWER_MASK 0x7fffffffUL /* least significant r bits */

/* initializes mt[N] with a seed */
void mt_ctx_init(mt_ctx *ctx, unsigned long s)
{
    int mti;

    ctx->mt[0] = s & 0xffffffffUL;
    for (mti = 1; mti < N; mti++) {
        ctx->mt[mti] =
            (1812433253UL * (ctx->mt[mti - 1] ^ (ctx->mt[mti - 1] >> 30))
             + mti);
        /* See Knuth TAOCP Vol2. 3rd Ed. P.106 for multiplier. */
        /* In the previous versions, MSBs of the seed affect   */
        /* only MSBs of the array mt[].                        */
        /* 2002/01/09 modified by Makoto Matsumoto             */
        ctx->mt[mti] &= 0xffffffffUL;
        /* for >32 bit machines */
    }
    ctx->mti = mti;
}

/* generates a random number on [0,0xffffffff]-interval */
unsigned long mt_ctx_int32(mt_ctx *ctx)
{
    unsigned long y;
    static unsigned long mag01[2] = {0x0UL, MATRIX_A};
    /* mag01[x] = x * MATRIX_A  for x=0,1 */

    if (ctx->mti >= N) { /* generate N words at one time */
        int kk;

        for (kk = 0; kk < N - M; kk++) {
            y = (ctx->mt[kk] & UPPER_MASK) | (ctx->mt[kk + 1] & LOWER_MASK);
            ctx->mt[kk] = ctx->mt[kk + M] ^ (y >> 1) ^ mag01[y & 0x1UL];
        }
        for (; kk < N - 1; kk++) {
            y = (ctx->mt[kk] & UPPER_MASK) | (ctx->mt[kk + 1] & LOWER_MASK);
            ctx->mt[kk] = ctx->mt[kk + (M - N)] ^ (y >> 1) ^ mag01[y & 0x1UL];
        }
        y = (ctx->mt[N - 1] & UPPER_MASK) | (ctx->mt[0] & LOWER_MASK);
        ctx->mt[N - 1] = ctx->mt[M - 1] ^ (y >> 1) ^ mag01[y & 0x1UL];

        ctx->mti = 0;
    }

    y = ctx->mt[ctx->mti++];

    /* Tempering */
    y ^= (y >> 11);
    y ^= (y << 7) & 0x9d2c5680UL;
    y ^= (y << 15) & 0xefc60000UL;
    y ^= (y >> 18);

    return y;
}

/* generates a random number on [0,1)-real-interval */
double mt_ctx_real2(mt_ctx *ctx)
{
    return mt_ctx_int32(ctx) * (1.0 / 4294967296.0);
    /* divided by 2^32 */
}
