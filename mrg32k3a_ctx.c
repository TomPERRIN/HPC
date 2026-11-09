/* ------------------------------------------------------------------- */
/* mrg32k3a_ctx.c                                                      */
/*                                                                     */
/* Thread-safe (context based) version of MRG32k3a of Pierre L'Ecuyer. */
/* Parameters come from MRG.h ; the algorithm is unchanged from MRG.c, */
/* only the state is a context given as a parameter.                   */
/* ------------------------------------------------------------------- */

#include "MRG.h"
#include "mrg32k3a_ctx.h"

void mrg32k3a_ctx_init(mrg32k3a_ctx *ctx, unsigned long long seed)
{
    unsigned long long s;

    if (seed == 0) seed = 1;          /* forbid the zero state */

    s = seed % (unsigned long long)MRG32K3A_M1;
    ctx->s10 = (s == 0) ? 1.0 : (double)s;
    ctx->s11 = 0.0;   ctx->s12 = 0.0;

    s = seed % (unsigned long long)MRG32K3A_M2;
    ctx->s20 = (s == 0) ? 1.0 : (double)s;
    ctx->s21 = 0.0;   ctx->s22 = 0.0;
}

double mrg32k3a_ctx_uniform(mrg32k3a_ctx *ctx)
{
    long   k;
    double p1, p2;
    /* Component 1 */
    p1 = MRG32K3A_A12 * ctx->s11 - MRG32K3A_A13N * ctx->s10;
    k = p1 / MRG32K3A_M1;   p1 -= k * MRG32K3A_M1;
    if (p1 < 0.0) p1 += MRG32K3A_M1;
    ctx->s10 = ctx->s11;   ctx->s11 = ctx->s12;   ctx->s12 = p1;
    /* Component 2 */
    p2 = MRG32K3A_A21 * ctx->s22 - MRG32K3A_A23N * ctx->s20;
    k  = p2 / MRG32K3A_M2;  p2 -= k * MRG32K3A_M2;
    if (p2 < 0.0) p2 += MRG32K3A_M2;
    ctx->s20 = ctx->s21;   ctx->s21 = ctx->s22;   ctx->s22 = p2;
    /* Combination */
    if (p1 <= p2) return ((p1 - p2 + MRG32K3A_M1) * MRG32K3A_NORM);
    else          return ((p1 - p2) * MRG32K3A_NORM);
}
