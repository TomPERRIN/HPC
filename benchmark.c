/* ------------------------------------------------------------------- */
/* benchmark.c                                                         */
/*                                                                     */
/* Performance test of the pseudorandom number generators:             */
/*   - rand()      the weak system generator (banned for science)      */
/*   - MT19937     Mersenne Twister, original code by Matsumoto        */
/*   - MRG32k3a    L'Ecuyer (MRG.h / MRG.c)                            */
/*   - MRG32k5a / MRG63k3a   high-quality but slower (bonus)           */
/*                                                                     */
/* Each generator is called N times (default 1 000 000 000 = 10^9)     */
/* drawing pseudo-random numbers between 0 and 1. A dummy sum (sink)   */
/* prevents the compiler from removing the loop.                       */
/*                                                                     */
/* Compile WITHOUT optimization first, then WITH it:                   */
/*   gcc benchmark.c MRG.c mt19937ar.c -o bench_O0                     */
/*   gcc -O2 benchmark.c MRG.c mt19937ar.c -o bench_O2                 */
/*                                                                     */
/* Set the number of calls and launch:                                 */
/*   ./bench_O0            (default: 10^9 calls)                       */
/*   ./bench_O2 1000000000 (same as default)                           */
/* ------------------------------------------------------------------- */

#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include "MRG.h"
#include "mt19937ar.h"

#define DEFAULT_CALLS 1000000000ULL

/** ------------------------------------------------------- *
  * now_sec   Monotonic wall-clock time in seconds
  *           (immune to system clock adjustments)
  * @return   time in seconds since an arbitrary origin
  * ------------------------------------------------------- **/

static double now_sec(void)
{
    struct timespec ts;

    clock_gettime(CLOCK_MONOTONIC, &ts);

    return (double)ts.tv_sec + 1.0e-9 * (double)ts.tv_nsec;
}

/** ------------------------------------------------------- *
  * report    Print one line of the results table
  *
  * @brief    Converts the elapsed time into the two usual
  *           performance metrics of a generator
  *
  * @param    name   generator name
  * @param    n      number of calls performed
  * @param    t      elapsed time in seconds
  * @param    sink   dummy sum, printed to show the loop
  *                  was not optimized away
  * ------------------------------------------------------- **/

static void report(const char *name, unsigned long long n,
                   double t, double sink)
{
    printf("%-12s %9.2f s   %8.3f ns/call   %10.1f Mcalls/s   (sink=%.4e)\n",
           name, t, t * 1.0e9 / (double)n, (double)n / t / 1.0e6, sink);
}

/** ------------------------------------------------------- *
  * main      Time N calls of each generator in [0..1)
  * ------------------------------------------------------- **/

int main(int argc, char *argv[])
{
    unsigned long long n, i;
    double t0, t, sink;
    const double RAND_NORM = 1.0 / (double)RAND_MAX;   /* avoid a division per call */

    n = DEFAULT_CALLS;
    if (argc > 1) n = strtoull(argv[1], NULL, 10);

    printf("Benchmark of pseudorandom number generators in [0..1)\n");
    printf("Calls per generator: %llu\n\n", n);

    /* --- 1) rand() : the weak system generator (reference) ---------- */
    t0 = now_sec();
    sink = 0.0;
    for (i = 0; i < n; i++)
        sink += (double)rand() * RAND_NORM;
    t = now_sec() - t0;
    report("rand()", n, t, sink);

    /* --- 2) MT19937 : Mersenne Twister (Matsumoto) ------------------- */
    init_genrand(5489UL);          /* original initialization routine */
    t0 = now_sec();
    sink = 0.0;
    for (i = 0; i < n; i++)
        sink += genrand_real2();   /* real version in [0..1) */
    t = now_sec() - t0;
    report("MT19937", n, t, sink);

    /* --- 3) MRG32k3a : L'Ecuyer -------------------------------------- */
    MRG32k3a_init(12345);
    t0 = now_sec();
    sink = 0.0;
    for (i = 0; i < n; i++)
        sink += MRG32k3a();
    t = now_sec() - t0;
    report("MRG32k3a", n, t, sink);

    /* --- 4) MRG32k5a : high quality but slower (bonus) ---------------- */
    MRG32k5a_init(12345);
    t0 = now_sec();
    sink = 0.0;
    for (i = 0; i < n; i++)
        sink += MRG32k5a();
    t = now_sec() - t0;
    report("MRG32k5a", n, t, sink);

    /* --- 5) MRG63k3a : 63-bit version (bonus) ------------------------- */
    MRG63k3a_init(12345);
    t0 = now_sec();
    sink = 0.0;
    for (i = 0; i < n; i++)
        sink += MRG63k3a();
    t = now_sec() - t0;
    report("MRG63k3a", n, t, sink);

    return 0;
}
