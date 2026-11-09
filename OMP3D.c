/* ------------------------------------------------------------------- */
/* OMP3D.c                                                             */
/*                                                                     */
/* Monte Carlo estimation of the volume of a sphere in 3D (4/3 PI R^3) */
/* with 95% confidence intervals around the simulated mean, computed   */
/* exactly as described in the appendix of the lab (Student law).      */
/*                                                                     */
/* The approach follows the 2D example OMP.c: the replicates are       */
/* computed in parallel with OpenMP, each replicate using its own      */
/* independent pseudorandom stream (thread-safe context based          */
/* generators, see mt19937ar_ctx.h and mrg32k3a_ctx.h).                */
/*                                                                     */
/* The user inputs:                                                    */
/*   - the number of replicates n (each replicate is an independent    */
/*     experiment, run under the same conditions with an independent   */
/*     random stream)                                                  */
/*   - the number of random drawings per replicate (sampling points)   */
/*   - the generator: MT19937 (Matsumoto) or MRG32k3a (L'Ecuyer)       */
/*                                                                     */
/* For each replicate: the points are drawn uniformly in the unit cube */
/* [0,1]^3, the fraction f falling inside the unit sphere (centered at */
/* the origin, radius R = 1) gives the volume estimate                 */
/*         Vi = 8 * f          (4/3 PI R^3, R = 1)                     */
/* and the PI estimate  PII = 6 * f  = 3/4 * Vi.                       */
/*                                                                     */
/* Statistics on the Xi = Vi (identical independent Gaussians):        */
/*   Xbar(n)  = sum(Xi)/n                       (arithmetical mean)    */
/*   S^2(n)   = sum((Xi-Xbar)^2)/(n-1)          (unbiased variance)    */
/*   SE       = S(n)/sqrt(n)  (standard error of the sample mean)      */
/*   R        = t(n-1, 0.975) * S(n)/sqrt(n)    (confidence radius)    */
/*   CI       = [Xbar - R, Xbar + R]            (95% confidence)       */
/*                                                                     */
/* All computations use double precision (only ~7 significant digits   */
/* in a float).                                                        */
/*                                                                     */
/* Compile with openMP flag:                                           */
/*   gcc -O2 -fopenmp OMP3D.c mrg32k3a_ctx.c mt19937ar_ctx.c -o mc3d   */
/*                                                                     */
/* Set the #of threads and launch: OMP_NUM_THREADS=8 ./mc3d            */
/* ------------------------------------------------------------------- */

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <omp.h>
#include "MRG.h"
#include "mrg32k3a_ctx.h"
#include "mt19937ar_ctx.h"

#define DEFAULT_POINTS 1000000000LL   /* 10^9 points per replicate */

enum { GEN_MT = 1, GEN_MRG = 2 };

/** ------------------------------------------------------- *
  * t975         Student law, t(n-1, 1-0.05/2) values
  *              (Table 1 of the appendix of the lab)
  *
  * @brief       Gives the confidence coefficient t0,975 for
  *              n experiments, used in R = t * S(n)/sqrt(n)
  *
  * @param  n    number of replicates (experiments)
  * @return      t(n-1, 0.975)
  * ------------------------------------------------------- **/

static double t975(int n)
{
    static const double table[30] = {
        12.706, 4.303, 3.182, 2.776, 2.571,   /*  1 -  5 */
         2.447, 2.365, 2.308, 2.262, 2.228,   /*  6 - 10 */
         2.201, 2.179, 2.160, 2.145, 2.131,   /* 11 - 15 */
         2.120, 2.110, 2.101, 2.093, 2.086,   /* 16 - 20 */
         2.080, 2.074, 2.069, 2.064, 2.060,   /* 21 - 25 */
         2.056, 2.052, 2.048, 2.045, 2.042    /* 26 - 30 */
    };

    if (n <= 30)   return table[n - 1];
    if (n <= 40)   return 2.021;
    if (n <= 80)   return 2.000;
    if (n <= 120)  return 1.980;
    return 1.960;                      /* n -> +oo */
}

/** ------------------------------------------------------- *
  * volume_estimate_mt    One replicate with MT19937
  *
  * @brief  Draws "points" triples (x,y,z) uniformly in [0,1]^3
  *         with an independent MT stream and returns the
  *         estimated sphere volume  Vi = 8 * in/points
  *
  * @param  points  number of sampling points (drawings/3)
  * @param  seed    seed of this replicate's stream
  * @return         the volume estimate Vi
  * ------------------------------------------------------- **/

static double volume_estimate_mt(long long points, unsigned long seed)
{
    mt_ctx ctx;
    long long in = 0, i;

    mt_ctx_init(&ctx, seed);

    for (i = 0; i < points; i++)
    {
        double x = mt_ctx_real2(&ctx);
        double y = mt_ctx_real2(&ctx);
        double z = mt_ctx_real2(&ctx);

        if (x * x + y * y + z * z <= 1.0)
            in++;
    }

    return 8.0 * (double)in / (double)points;
}

/** ------------------------------------------------------- *
  * volume_estimate_mrg   One replicate with MRG32k3a
  *                       (same as volume_estimate_mt)
  * ------------------------------------------------------- **/

static double volume_estimate_mrg(long long points,
                                  unsigned long long seed)
{
    mrg32k3a_ctx ctx;
    long long in = 0, i;

    mrg32k3a_ctx_init(&ctx, seed);

    for (i = 0; i < points; i++)
    {
        double x = mrg32k3a_ctx_uniform(&ctx);
        double y = mrg32k3a_ctx_uniform(&ctx);
        double z = mrg32k3a_ctx_uniform(&ctx);

        if (x * x + y * y + z * z <= 1.0)
            in++;
    }

    return 8.0 * (double)in / (double)points;
}

/** ------------------------------------------------------- *
  * main    Monte Carlo 3D + 95% confidence interval
  * ------------------------------------------------------- **/

int main(void)
{
    int      nRep, gen, r;
    long long points;
    double  *V;                 /* one volume estimate per replicate */
    double   mean, S2, S, SE, radius, tcoef;
    double   piMean, piRadius, exactV, exactPI;
    double   t0, t;

    printf("Monte Carlo 3D: volume of a sphere (4/3 PI R^3, R = 1)\n");

    /* --- user inputs ------------------------------------------------ */
    do {
        printf("Number of replicates n (>1): ");
        if (scanf("%d", &nRep) != 1) return 1;
    } while (nRep < 2);

    printf("Number of points per replicate [0 = default %lld]: ",
           DEFAULT_POINTS);
    if (scanf("%lld", &points) != 1) return 1;
    if (points <= 0) points = DEFAULT_POINTS;

    do {
        printf("Generator (1 = MT19937, 2 = MRG32k3a): ");
        if (scanf("%d", &gen) != 1) return 1;
    } while (gen != GEN_MT && gen != GEN_MRG);

    /* --- replicates in parallel, one independent stream each -------- */
    V = (double *)malloc((size_t)nRep * sizeof(double));
    if (V == NULL) return 1;

    printf("\nComputing %d replicates x %lld points with %d threads...\n",
           nRep, points, omp_get_max_threads());

    t0 = omp_get_wtime();

    #pragma omp parallel for
    for (r = 0; r < nRep; r++)
    {
        if (gen == GEN_MT)
            /* one MT stream per replicate, seeds far apart */
            V[r] = volume_estimate_mt(points, 5489UL + 10007UL * (unsigned long)r);
        else
            /* one MRG32k3a stream per replicate, seeds far apart */
            V[r] = volume_estimate_mrg(points, 12345ULL + 10007ULL * (unsigned long long)r);
    }

    t = omp_get_wtime() - t0;

    /* --- statistics on the Xi = Vi (appendix of the lab) ------------- */
    mean = 0.0;
    for (r = 0; r < nRep; r++)
        mean += V[r];
    mean /= (double)nRep;

    S2 = 0.0;
    for (r = 0; r < nRep; r++)
        S2 += (V[r] - mean) * (V[r] - mean);
    S2 /= (double)(nRep - 1);           /* unbiased estimate of sigma^2 */

    S      = sqrt(S2);
    SE     = S / sqrt((double)nRep);    /* standard error of the mean */
    tcoef  = t975(nRep);
    radius = tcoef * SE;                /* 95% confidence radius */

    piMean   = 0.75 * mean;             /* PI = 3/4 * volume (R = 1) */
    piRadius = 0.75 * radius;
    exactV   = 4.0 / 3.0 * acos(-1.0);  /* 4/3 PI, R = 1 */
    exactPI  = acos(-1.0);

    /* --- report ------------------------------------------------------ */
    printf("\nPer-replicate volume estimates\n");
    for (r = 0; r < nRep; r++)
        printf("Replicate %2d: V = %.10f   PI = %.10f\n",
               r + 1, V[r], 0.75 * V[r]);

    printf("\n----- Statistics (n = %d) -----\n", nRep);
    printf("Mean volume   Xbar(n):    %.10f\n", mean);
    printf("Variance      S^2(n):     %.10e\n", S2);
    printf("Std deviation S(n):       %.10e\n", S);
    printf("Standard error S(n)/sqrt(n): %.10e\n", SE);
    printf("Student       t(n-1,0.975): %.3f\n", tcoef);
    printf("Confidence radius R:      %.10e\n\n", radius);

    printf("Exact volume (4/3 PI):    %.10f\n", exactV);
    printf("Estimated volume:         %.10f\n", mean);
    printf("95%% CI on volume:         [%.10f , %.10f]\n",
           mean - radius, mean + radius);
    printf("Exact value inside CI:    %s\n\n",
           (mean - radius <= exactV && exactV <= mean + radius) ? "yes" : "no");

    printf("Exact PI:                 %.10f\n", exactPI);
    printf("Estimated PI:             %.10f\n", piMean);
    printf("95%% CI on PI:             [%.10f , %.10f]\n",
           piMean - piRadius, piMean + piRadius);
    printf("Exact value inside CI:    %s\n\n",
           (piMean - piRadius <= exactPI && exactPI <= piMean + piRadius) ? "yes" : "no");

    printf("Absolute error (volume):  %.10e\n", fabs(mean - exactV));
    printf("Elapsed time:             %.2f s\n", t);

    free(V);

    return 0;
}
