
#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <math.h>

#define N 10

double simulPi(int nbPoints) {
    int i;
    double xr = 0.;
    double yr = 0.;
    double pi = 0.;
    int inDisk = 0;

    for (i = 0; i < nbPoints; i++) {
        xr = (double)rand() / RAND_MAX;
        yr = (double)rand() / RAND_MAX;

        if ((xr * xr + yr * yr) <= 1) {
            inDisk++;
        }
    }

    pi = ((double)inDisk / (double)nbPoints) * 4.;

    return pi;
}

void experiment(int nbPoints) {
    double results[N];
    double meanPI = 0.;
    double absoluteError;
    double relativeError;
    int i;

    for (i = 0; i < N; i++) {
        results[i] = simulPi(nbPoints);
        meanPI += results[i];
    }

    meanPI /= N;

    absoluteError = fabs(meanPI - M_PI);
    relativeError = absoluteError / M_PI;

    printf("\nNombre de points : %d\n", nbPoints);

    for (i = 0; i < N; i++) {
        printf("Experience %d : %.12f\n", i + 1, results[i]);
    }

    printf("Moyenne de PI       : %.12f\n", meanPI);
    printf("M_PI                : %.12f\n", M_PI);
    printf("Erreur absolue      : %.12f\n", absoluteError);
    printf("Erreur relative     : %.12f\n", relativeError);
}

int main() {

    srand((unsigned int)time(NULL));

    experiment(1000);
    experiment(1000000);
    experiment(1000000000);

    return 0;
}
