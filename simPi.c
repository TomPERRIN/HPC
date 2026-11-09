
#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#define NB_POINTS 1000000000

double simulPi(int nbPoints) {
    int i;
    double xr, yr;
    int inDisk = 0;

    for (i = 0; i < nbPoints; i++) {
        xr = (double)rand() / RAND_MAX;
        yr = (double)rand() / RAND_MAX;

        if ((xr * xr + yr * yr) <= 1.0) {
            inDisk++;
        }
    }

    return 4.0 * (double)inDisk / (double)nbPoints;
}

int main(int argc, char *argv[]) {

    // Vérification de l'argument
    if (argc != 2) {
        fprintf(stderr, "Usage : %s seed\n", argv[0]);
        return 1;
    }

    // Récupération de la graine
    unsigned int seed = (unsigned int)atoi(argv[1]);

    // Initialisation du générateur aléatoire
    srand(seed);

    // Simulation de PI
    double pi = simulPi(NB_POINTS);

    // Affichage sur la sortie standard
    printf("%.12f\n", pi);

    return 0;
}
