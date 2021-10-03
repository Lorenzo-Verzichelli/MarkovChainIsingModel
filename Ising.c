#include <stdio.h>
#include <math.h>
#include <stdlib.h>

#include "Ran2.c"

#define Nlato 50
#define h 0.1
#define beta 0.1
#define Nmisure 1000000


//calcola la magnetizzazione: field è la matrice degli spin//
int magnetiz(int field[Nlato][Nlato]) {
    int m = 0;
    for (int i=0; i<Nlato; i++){
        for (int j=0; j<Nlato; j++){
            m += field[i][j];
        }
    }
    return m;
}

//inizializza con tutti gli spin su, che segniamo con 1//
void init_freddo(int field[Nlato][Nlato]){
    for (int i=0; i<Nlato; i++){
        for (int j=0; j<Nlato; j++){
            field[i][j] = 1;
        }
    }
}

//inizializzi a caldo, ovvero con ogni spin a caso//
void init_casuale(int field[Nlato][Nlato], long* seed) {
    for (int i=0; i<Nlato; i++){
        for (int j=0; j<Nlato; j++){
            field[i][j] = 1 - 2*((int)(2*ran2(seed)));
        }
    }
}

//calcolo la forza sul sito i j del reticolo
int forza(int field[Nlato][Nlato], int i, int j){
	return field[(i+1)%Nlato][j] + field[(i+Nlato-1)%Nlato][j] + field[i][(j+1)%Nlato] + field[i][(j+Nlato-1)%Nlato];
}

//faccio un passo del metropolis
void passo_metropolis(int field[Nlato][Nlato], long* seed) {
    int i, j, f;
    i = (int)(Nlato*ran2(seed));
    j = (int)(Nlato*ran2(seed));
    f = forza(field, i, j);
    if (ran2(seed) < exp(-2*beta*field[i][j]*(f+h))) field[i][j] = -field[i][j];
}

//calcolo l'energia totale
double energia(int field[Nlato][Nlato]) {
    int u = 0;
    for (int i=0; i<Nlato; i++){
        for (int j=0; j<Nlato; j++){
            u += field[i][j]*forza(field, i, j);
            //printf("u = %d\n", u);
        }
    }
    return -0.5*u - magnetiz(field)*h;
}

//stampa la matrice di spin nella configurazione corrente sullo schermo//
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
    init_casuale(field, &seed);
    //print_field(field);
    printf("%d \n", magnetiz(field));
    for (int i=0; i<Nmisure; i++){
        passo_metropolis(field, &seed);
        //print_field(field);
    }
    //print_field(field);
    printf("magnetiz: %d \n", magnetiz(field));
    printf("energia : %f \n", energia(field));
    printf("%d\n", forza(field, 0, 0));
    return 0;
}



