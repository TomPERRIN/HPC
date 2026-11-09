/* ------------------------------------------------------------------- */
/* mt19937ar_ctx.h                                                     */
/*                                                                     */
/* Thread-safe (context based) version of MT19937.                     */
/*                                                                     */
/* The generation/initialization algorithms are UNCHANGED from the     */
/* original mt19937ar.c/.h of Takuji Nishimura and Makoto Matsumoto    */
/* (BSD license, see mt19937ar.c). Only difference: the state is no    */
/* longer a global variable but a context given as a parameter, so     */
/* that each OpenMP thread / each replicate can own an independent     */
/* stream (like xorshift64next(&state) in OMP.c).                      */
/*                                                                     */
/* Compile with:  gcc -c mt19937ar_ctx.c                               */
/* ------------------------------------------------------------------- */

#ifndef MT19937AR_CTX_H
#define MT19937AR_CTX_H

#define MT_CTX_N 624

typedef struct
{
    unsigned long mt[MT_CTX_N];   /* the array for the state vector */
    int           mti;            /* mti==MT_CTX_N+1 means uninitialized */
} mt_ctx;

/* initializes the context with a seed (same algorithm as init_genrand) */
void mt_ctx_init(mt_ctx *ctx, unsigned long seed);

/* generates a random unsigned 32-bit integer (same as genrand_int32) */
unsigned long mt_ctx_int32(mt_ctx *ctx);

/* generates a random number on [0,1)-real-interval (same as genrand_real2) */
double mt_ctx_real2(mt_ctx *ctx);

#endif /* MT19937AR_CTX_H */
