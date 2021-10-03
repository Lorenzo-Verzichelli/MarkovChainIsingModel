#include <stdio.h>
#include <math.h>
#include <stdlib.h>

#include "Ran2.c"

#define Nlato 10
#define h 0.1
#define beta 0.1
#define Nmisure 1000000

int magnetiz(int field[Nlato][Nlato]) {
    int m = 0;
    for (int i=0; i<Nlato; i++){
        for (int j=0; j<Nlato; j++){
            m += field[i][j];
        }
    }
    return m;
}

void init_freddo(int field[Nlato][Nlato]){
    for (int i=0; i<Nlato; i++){
        for (int j=0; j<Nlato; j++){
            field[i][j] = 1;
        }
    }
}

void init_casuale(int field[Nlato][Nlato], long* seed) {
    for (int i=0; i<Nlato; i++){
        for (int j=0; j<Nlato; j++){
            field[i][j] = 1 - 2*((int)(2*ran2(seed)));
        }
    }
}

void passo_metropolis(int field[Nlato][Nlato], long* seed) {
    int i, j, f;
    i = (int)(Nlato*ran2(seed));
    j = (int)(Nlato*ran2(seed));
    f = field[(i+1)%Nlato][j] + field[(i-1)%Nlato][j]+field[i][(j+1)%Nlato]+field[i][(j-1)%Nlato];
    if (ran2(seed) < exp(-2*beta*field[i][j]*(f+h))) field[i][j] = -field[i][j];
}

void print_field(int field[Nlato][Nlato]) {
    for (int i=0; i<Nlato; i++){
        for (int j= 0; j<Nlato; j++) {
            if (field[i][j]+1) printf(" 1 ");
            else printf("-1 ");
        }
        printf("\n");
    }
    printf("\n");
}

int main(void) {
    long seed=-3;
    int field[Nlato][Nlato];
    init_freddo(field);
    print_field(field);
    printf("%d \n", magnetiz(field));
    for (int i=0; i<Nmisure; i++){
        passo_metropolis(field, &seed);
        //print_field(field);
    }
    print_field(field);
    printf("%d \n", magnetiz(field));
    return 0;
}
