#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#define MAX_CORR 5000

void statistica(double temp, int n, double* med, double* M2, double* M3, double* M4) {
    double delta = temp - *med;
    double delta_n = delta/n;
    *med = *med + delta_n;
    double term1 = delta*delta_n*(n - 1);
    if (M3 != NULL) {
        double delta_n2 = delta_n * delta_n;
        *M4 = *M4 + term1*delta_n2*(n*n - 3*n + 3) + 6*delta_n2*(*M2) - 4*delta_n*(*M3);
        *M3 = *M3 + term1*delta_n*(n - 2) - 3*delta_n*(*M2);
    }
    *M2 = *M2 + term1;
}

int main (int argc, char *argv[]) {
    FILE *input1, *output;
    switch (argc) {
        case 1: input1 = fopen("Default_out.txt", "r");
                output = fopen("Default_cor.txt", "w"); break;
        case 2: input1 = fopen(argv[1], "r");
                output = fopen("Default_cor.txt", "w"); break;
        default: input1 = fopen(argv[1], "r");
                 output = fopen(argv[2], "w");
                 perror("ignorati tutti gli argomenti, tranne i primi 2");
    }
    if (input1 == NULL) {
        perror("fopen input1");
        return EXIT_FAILURE;
    }

    char c = fgetc(input1);
    if (c == '#') {
        while (fgetc(input1) != '\n');
    }
    else ungetc(c, input1);
    fpos_t inizio; fgetpos(input1, &inizio);

    int tempE, tempM, cont = 0;
    double medE=0, medM=0, M2E=0, M2M=0, M3M=0, M4M=0;
    while (fscanf(input1, "%d\t%d\n", &tempE, &tempM) != EOF){
        cont ++;
        statistica((double)tempE, cont, &medE, &M2E, NULL, NULL);
        statistica((double)tempM, cont, &medM, &M2M, &M3M, &M4M);
    }
    printf("energia: %f\nmagnetizzazione: %f\n", medE, medM);
    printf("%d, %e, %e, %e\n", cont, M2E, M2M, M4M);

    int *misE, *misM;
    misE = malloc(cont * sizeof(int));
    misM = malloc(cont * sizeof(int));

    fsetpos(input1, &inizio);
    printf("ciao!");
    for (int i = 0; i<cont; i++) {
        fscanf(input1, "%d\t%d\n", &tempE, &tempM);
        misE[i] = tempE - medE;
        misM[i] = tempM - medM;
    }

    double corrE, sumE = 0, corrM, sumM = 0, varE = M2E/cont, varM = M2M/cont;
    for (int k = 0; k <= MAX_CORR; k++) {
        corrE = 0;
        corrM = 0;
        for (int i = 0; i < cont-k; i++) {
            corrE += misE[i]*misE[i+k];
            corrM += misM[i]*misM[i+k];
        }
        corrE /= varE*(cont - k);
        corrM /= varM*(cont - k);
        sumE += corrE;
        sumM += corrM;
        fprintf(output, "%f\t%f\t%f\t%f\n", corrE, sumE, corrM, sumM);
    }
    fclose(input1); fclose(output);
    return 0;
}
