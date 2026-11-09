/* ------------------------------------------------------------------- */
/* openMPcomputePI                                                     */
/*                                                                     */
/* Compile with openMP flag: gcc -O2 -fopenmp monte_carlo_pi.c -o mcPI */
/*                                                                     */
/* Set the #of threads and launch: OMP_NUM_THREADS=4 ./mcPI            */
/* ------------------------------------------------------------------- */

#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <math.h>
#include <omp.h>

#define NUM_REPLICATES 8
#define NUM_POINTS 10000000

/** ------------------------------------------------------- *
  * xorshift64next       Simple thread-local pseudorandom   *
  *                      number generator                   *
  *                                                         *
  * @brief  This function an unsigned int random generator  *
  *         it is thread safe thanks to its state parameter *
  *                                                         *
  * @param  inState   The generator state                   *
  *                                                         *
  * @return a pseudorandom unsigned int                     *
  * ------------------------------------------------------- **/

static inline uint64_t xorshift64next(uint64_t *inState)
{
    uint64_t x = *inState;
	
    x ^= x >> 12;
    x ^= x << 25;
    x ^= x >> 27;
	
    *inState = x;
    
	return x * UINT64_C(2685821657736338717);
}

/** ------------------------------------------------------- *
  * uniform01      Give a pseudorandom number in [0..1 ]    *
  *                                                         *
  * @brief  This function takes an "random" integer         *
  *         coming form xorshift64next it supresses         *
  *         exessive bits (keeps enoough for a double       *
  *         mantissa up to 53 bits                          *
  *         (1 sign + 10 bits exponent + 53 for mantissa    *   
  *                                                         *
  * @param  inState   The generator state                   *
  *                                                         *
  * @return a pseudorandom number between 0 & 1             *
  * ------------------------------------------------------- **/

static inline double uniform01(uint64_t * inState)
{
    return (xorshift64next(inState) >> 11) * 0x1.0p-53;
}

/** ------------------------------------------------------- *
  * main 	Gives an estimate of PI                         *
  *         Computes the replicates in parallel             *
  *         with a number of OpenMP threads                 *
  * ------------------------------------------------------- **/

int main(void)
{
    double estimates[NUM_REPLICATES];

    #pragma omp parallel for
	for (int r = 0; r < NUM_REPLICATES; r++) 
	{
        // the state of the generator depends 
		// on the replicate number
	
        uint64_t state =
            UINT64_C(123456789) +
            (uint64_t)(r + 1) * UINT64_C(987654321);

        long long inDisk = 0;

        for (long long i = 0; i < NUM_POINTS; i++) 
		{
            double x = uniform01(&state);
            double y = uniform01(&state);

            if (x * x + y * y <= 1.0) 
			{
               inDisk++;
            }
        }

        estimates[r] = 4.0 * (double)inDisk / NUM_POINTS;
    }

    /* Combine the independent replicate estimates */
    double mean = 0.0;

    printf("Monte Carlo estimation of PI\n");
    printf("Replicates: %d\n", NUM_REPLICATES);
    printf("Number of points per replicate: %d\n\n", NUM_POINTS);

    for (int r = 0; r < NUM_REPLICATES; r++) 
	{
        printf("Replicate %2d: %.10f\n", r + 1, estimates[r]);
    
		mean += estimates[r];
    }

    mean /= NUM_REPLICATES;

    printf("\nMean estimate: %.10f\n", mean);
    printf("Exact PI:      %.10f\n", acos(-1.0));
    printf("Absolute error: %.10e\n", fabs(mean - acos(-1.0)));

    return 0;
}
