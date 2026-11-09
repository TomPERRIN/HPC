/* ------------------------------------------------------------------- */
/* MRG.h                                                               */
/*                                                                     */
/* Combined Multiple Recursive Generators of Pierre L'Ecuyer           */
/*   - MRG32k3a : 32-bit components, period ~ 2^191                    */
/*   - MRG32k5a : 32-bit components, period ~ 2^319                    */
/*   - MRG63k3a : 63-bit components, period ~ 2^377                    */
/*                                                                     */
/* Reference: P. L'Ecuyer, "Good Parameter Sets for Combined Multiple  */
/* Recursive Random Number Generators", Operations Research, 47(1),    */
/* 1999, pp 159-164.                                                   */
/*                                                                     */
/* Each generator must be initialised with its own ..._init() function */
/* before its first call: the all-zero state is a fixed point of the   */
/* generator and must never be used.                                   */
/*                                                                     */
/* Compile with:  gcc -c MRG.c                                         */
/* ------------------------------------------------------------------- */

#ifndef MRG_H
#define MRG_H

/* ------------------------- MRG32k3a -------------------------------- */
/* MRG(3, [0, a12, -a13n], m1 ; [a21, 0, -a23n], m2)                   */
/* norm = 1 / (m1 + 1)                                                 */
#define MRG32K3A_NORM  2.328306549295728e-10
#define MRG32K3A_M1    4294967087.0
#define MRG32K3A_M2    4294944443.0
#define MRG32K3A_A12      1403580.0
#define MRG32K3A_A13N      810728.0
#define MRG32K3A_A21       527612.0
#define MRG32K3A_A23N     1370589.0

/* ------------------------- MRG32k5a -------------------------------- */
/* MRG(5, [0, 0, a12, 0, -a15n], m1 ; [a21, 0, a23, 0, -a25n], m2)     */
/* norm = 1 / (m1 + 1)                                                 */
#define MRG32K5A_NORM  2.3283163396834613e-10
#define MRG32K5A_M1    4294949027.0
#define MRG32K5A_M2    4294934327.0
#define MRG32K5A_A12      1154721.0
#define MRG32K5A_A14     1739991.0
#define MRG32K5A_A15N    1108499.0
#define MRG32K5A_A21     1776413.0
#define MRG32K5A_A23      865203.0
#define MRG32K5A_A25N    1641052.0

/* ------------------------- MRG63k3a -------------------------------- */
/* MRG(3, [0, a12, -a13n], m1 ; [a21, 0, -a23n], m2) on 63 bits        */
/* norm = 1 / (m1 + 1)                                                 */
#define MRG63K3A_NORM  1.0842021724855052e-19
#define MRG63K3A_M1    9223372036854769163LL
#define MRG63K3A_M2    9223372036854754679LL
#define MRG63K3A_A12      1754669720
#define MRG63K3A_Q12      5256471877
#define MRG63K3A_R12       251304723
#define MRG63K3A_A13N     3182104042
#define MRG63K3A_Q13      2898513661
#define MRG63K3A_R13       394451401
#define MRG63K3A_A21     31387477935
#define MRG63K3A_Q21       293855150
#define MRG63K3A_R21       143639429
#define MRG63K3A_A23N     6199136374
#define MRG63K3A_Q23      1487847900
#define MRG63K3A_R23       985240079

/* ------------------------ Prototypes ------------------------------- */
void  MRG32k3a_init (unsigned long long seed);
double MRG32k3a     (void);

void  MRG32k5a_init (unsigned long long seed);
double MRG32k5a     (void);

void  MRG63k3a_init (unsigned long long seed);
double MRG63k3a     (void);

#endif /* MRG_H */
