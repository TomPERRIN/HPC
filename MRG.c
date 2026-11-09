/* ------------------------------------------------------------------- */
/* MRG.c                                                               */
/*                                                                     */
/* Combined Multiple Recursive Generators of Pierre L'Ecuyer           */
/* (MRG32k3a, MRG32k5a, MRG63k3a). The parameters (defines) are in     */
/* MRG.h ; each generator has its own state and its own                */
/* initialisation function (the all-zero state is forbidden).          */
/*                                                                     */
/* Compile with:  gcc -O2 -c MRG.c                                     */
/* ------------------------------------------------------------------- */

#include "MRG.h"

/* ------------------------------------------------------------------- */
/* MRG32k3a                                                            */
/* ------------------------------------------------------------------- */

static double s10_3a, s11_3a, s12_3a;
static double s20_3a, s21_3a, s22_3a;

void MRG32k3a_init (unsigned long long seed)
{
    unsigned long long s;

    if (seed == 0) seed = 1;          /* forbid the zero state */

    s = seed % (unsigned long long)MRG32K3A_M1;
    s10_3a = (s == 0) ? 1.0 : (double)s;
    s11_3a = 0.0;   s12_3a = 0.0;

    s = seed % (unsigned long long)MRG32K3A_M2;
    s20_3a = (s == 0) ? 1.0 : (double)s;
    s21_3a = 0.0;   s22_3a = 0.0;
}

double MRG32k3a (void)
{
    long   k;
    double p1, p2;
    /* Component 1 */
    p1 = MRG32K3A_A12 * s11_3a - MRG32K3A_A13N * s10_3a;
    k = p1 / MRG32K3A_M1;   p1 -= k * MRG32K3A_M1;
    if (p1 < 0.0) p1 += MRG32K3A_M1;
    s10_3a = s11_3a;   s11_3a = s12_3a;   s12_3a = p1;
    /* Component 2 */
    p2 = MRG32K3A_A21 * s22_3a - MRG32K3A_A23N * s20_3a;
    k  = p2 / MRG32K3A_M2;  p2 -= k * MRG32K3A_M2;
    if (p2 < 0.0) p2 += MRG32K3A_M2;
    s20_3a = s21_3a;   s21_3a = s22_3a;   s22_3a = p2;
    /* Combination */
    if (p1 <= p2) return ((p1 - p2 + MRG32K3A_M1) * MRG32K3A_NORM);
    else          return ((p1 - p2) * MRG32K3A_NORM);
}

/* ------------------------------------------------------------------- */
/* MRG32k5a                                                            */
/* ------------------------------------------------------------------- */

static double s10_5a, s11_5a, s12_5a, s13_5a, s14_5a;
static double s20_5a, s21_5a, s22_5a, s23_5a, s24_5a;

void MRG32k5a_init (unsigned long long seed)
{
    unsigned long long s;

    if (seed == 0) seed = 1;          /* forbid the zero state */

    s = seed % (unsigned long long)MRG32K5A_M1;
    s10_5a = (s == 0) ? 1.0 : (double)s;
    s11_5a = 0.0;   s12_5a = 0.0;   s13_5a = 0.0;   s14_5a = 0.0;

    s = seed % (unsigned long long)MRG32K5A_M2;
    s20_5a = (s == 0) ? 1.0 : (double)s;
    s21_5a = 0.0;   s22_5a = 0.0;   s23_5a = 0.0;   s24_5a = 0.0;
}

double MRG32k5a (void)
{
    long   k;
    double p1, p2;
    /* Component 1 */
    p1 = MRG32K5A_A12 * s13_5a - MRG32K5A_A15N * s10_5a;
    if (p1 > 0.0) p1 -= MRG32K5A_A14 * MRG32K5A_M1;
    p1 += MRG32K5A_A14 * s11_5a;
    k = p1 / MRG32K5A_M1;    p1 -= k * MRG32K5A_M1;
    if (p1 < 0.0) p1 += MRG32K5A_M1;
    s10_5a = s11_5a;   s11_5a = s12_5a;   s12_5a = s13_5a;
    s13_5a = s14_5a;   s14_5a = p1;
    /* Component 2 */
    p2 = MRG32K5A_A21 * s24_5a - MRG32K5A_A25N * s20_5a;
    if (p2 > 0.0) p2 -= MRG32K5A_A23 * MRG32K5A_M2;
    p2 += MRG32K5A_A23 * s22_5a;
    k  = p2 / MRG32K5A_M2;   p2 -= k * MRG32K5A_M2;
    if (p2 < 0.0) p2 += MRG32K5A_M2;
    s20_5a = s21_5a;   s21_5a = s22_5a;   s22_5a = s23_5a;
    s23_5a = s24_5a;   s24_5a = p2;
    /* Combination */
    if (p1 <= p2) return ((p1 - p2 + MRG32K5A_M1) * MRG32K5A_NORM);
    else          return ((p1 - p2) * MRG32K5A_NORM);
}

/* ------------------------------------------------------------------- */
/* MRG63k3a                                                            */
/* ------------------------------------------------------------------- */

static long long s10_63, s11_63, s12_63;
static long long s20_63, s21_63, s22_63;

void MRG63k3a_init (unsigned long long seed)
{
    if (seed == 0) seed = 1;          /* forbid the zero state */

    s10_63 = (long long)(seed % (unsigned long long)MRG63K3A_M1);
    if (s10_63 == 0) s10_63 = 1;
    s11_63 = 0;   s12_63 = 0;

    s20_63 = (long long)(seed % (unsigned long long)MRG63K3A_M2);
    if (s20_63 == 0) s20_63 = 1;
    s21_63 = 0;   s22_63 = 0;
}

double MRG63k3a (void)
{
    long long h, p12, p13, p21, p23;
    /* Component 1 */
    h = s10_63 / MRG63K3A_Q13;
    p13 = MRG63K3A_A13N * (s10_63 - h * MRG63K3A_Q13) - h * MRG63K3A_R13;
    h = s11_63 / MRG63K3A_Q12;
    p12 = MRG63K3A_A12  * (s11_63 - h * MRG63K3A_Q12) - h * MRG63K3A_R12;
    if (p13 < 0) p13 += MRG63K3A_M1;
    if (p12 < 0) p12 += MRG63K3A_M1 - p13;  else p12 -= p13;
    if (p12 < 0) p12 += MRG63K3A_M1;
    s10_63 = s11_63;  s11_63 = s12_63;  s12_63 = p12;
    /* Component 2 */
    h = s20_63 / MRG63K3A_Q23;
    p23 = MRG63K3A_A23N * (s20_63 - h * MRG63K3A_Q23) - h * MRG63K3A_R23;
    h = s22_63 / MRG63K3A_Q21;
    p21 = MRG63K3A_A21  * (s22_63 - h * MRG63K3A_Q21) - h * MRG63K3A_R21;
    if (p23 < 0) p23 += MRG63K3A_M2;
    if (p21 < 0) p21 += MRG63K3A_M2 - p23;  else p21 -= p23;
    if (p21 < 0) p21 += MRG63K3A_M2;
    s20_63 = s21_63;  s21_63 = s22_63;  s22_63 = p21;
    /* Combination */
    if (p12 > p21) return ((p12 - p21) * MRG63K3A_NORM);
    else           return ((p12 - p21 + MRG63K3A_M1) * MRG63K3A_NORM);
}
