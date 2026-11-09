/* ------------------------------------------------------------------- */
/* mrg32k3a_ctx.h                                                      */
/*                                                                     */
/* Thread-safe (context based) version of MRG32k3a of Pierre L'Ecuyer  */
/* (parameters defined in MRG.h, algorithm unchanged from MRG.c).      */
/*                                                                     */
/* Only difference with MRG.c: the state is no longer a set of static  */
/* variables but a context given as a parameter, so that each OpenMP   */
/* thread / each replicate can own an independent stream.              */
/*                                                                     */
/* Reference: P. L'Ecuyer, "Good Parameter Sets for Combined Multiple  */
/* Recursive Random Number Generators", Operations Research, 47(1),    */
/* 1999, pp 159-164.                                                   */
/*                                                                     */
/* Compile with:  gcc -c mrg32k3a_ctx.c                                */
/* ------------------------------------------------------------------- */

#ifndef MRG32K3A_CTX_H
#define MRG32K3A_CTX_H

typedef struct
{
    double s10, s11, s12;   /* component 1 */
    double s20, s21, s22;   /* component 2 */
} mrg32k3a_ctx;

/* initializes the context with a seed (the all-zero state is forbidden,
   as in MRG32k3a_init of MRG.c) */
void mrg32k3a_ctx_init(mrg32k3a_ctx *ctx, unsigned long long seed);

/* generates a random number on [0,1)-real-interval (same as MRG32k3a) */
double mrg32k3a_ctx_uniform(mrg32k3a_ctx *ctx);

#endif /* MRG32K3A_CTX_H */
