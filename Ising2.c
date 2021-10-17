#include <stdio.h>
#include <math.h>
#include <stdlib.h>

#include "Ran2.c"

#define Ngeo 5

void geometry(int geo[][Ngeo], int Nlato) {
    for (int i=0; i<Nlato-1; i++) {
        geo[i][1] = i+1;
        geo[i+1][0] = i;
        geo[i+1][2] = Nlato*i;
        geo[i][3] = Nlato*(i+1);
        geo[i][4] = Nlato*i;
    }
    geo[0][0] = Nlato-1;
    geo[Nlato-1][1] = 0;
    geo[0][2] = Nlato*(Nlato-1);
    geo[Nlato-1][3] = 0;
    geo[Nlato-1][4] = Nlato*(Nlato-1);
}
void init_file(int *field, int Nlato, FILE *init_f, long *seed) {
    fscanf(init_f, "%d\n", seed);
    for (int i=0; i<Nlato; i++){
        for (int j=0; j<Nlato-1; j++) {
            fscanf(init_f, "%d ", field+i*Nlato+j);
        }
        fscanf(init_f, "%d\n", field+i*Nlato+Nlato-1);
    }
}

void init_caldo(int *field, long *seed, int Nlato) {
    //inizializza il reticolo a caso
    for (int i=0; i<Nlato; i++){
        for (int j=0; j<Nlato; j++){
            field[i*Nlato + j] = 1 - 2*(ranInt(seed)%2);
        }
    }
}

void init_freddo(int *field, int Nlato) {
    //inizializza il reticolo con tutti gli spin in su
    for (int i=0; i<Nlato; i++){
        for (int j=0; j<Nlato; j++){
            field[i*Nlato + j] = 1;
        }
    }
}

void print_field(int *field, int Nlato, FILE *output){
    for (int i=0; i<Nlato; i++){
        for (int j= 0; j<Nlato; j++) {
            if (field[i*Nlato + j]+1) fprintf(output, " 1 ");
            else fprintf(output, "-1 ");
        }
        fprintf(output, "\n");
    }
    fprintf(output, "\n");
}

int forza(int *field, int Nlato, int i, int j, int geo[][Ngeo]){
	return field[geo[i][2] + j] + field[geo[i][3] + j] + field[geo[i][4] + geo[j][0]] + field[geo[i][4] + geo[j][1]];
}

void passo_metropolis(int *field, int Nlato, double beta, long* seed, int geo[][Ngeo]) {
    int i, j, f, s;
    long r = ranInt(seed);
    i = r % Nlato;
    j = (r/Nlato) % Nlato;
    s = field[geo[i][4]+j];
    f = forza(field, Nlato, i, j, geo);
    if (ran2(seed) < exp(-2*beta*s*(f))) field[geo[i][4]+j] = -s;
}

int energia(int *field, int Nlato, int geo[][Ngeo]) {
    int u = 0;
    for (int i=0; i<Nlato; i++){
        for (int j=0; j<Nlato; j++){
            u += field[geo[i][4] + j]*forza(field, Nlato, i, j, geo);
        }
    }
    return -u/2;
}

int tot_magnetizza(int *field, int Nlato) {
    int m = 0;
    for (int i=0; i<Nlato; i++){
        for (int j=0; j<Nlato; j++){
            m += field[i*Nlato + j];
        }
    }
    return m;
}

int main(int argc, char *argv[]) {
    //stampa help
    if(argc == 2 && argv[1][0] == '-' && argv[1][1] == 'h' && argv[1][2] == '\0') {
        printf("ARGOMENTI:\n - output file\n - beta\n - lato reticolo\n - opzione init (0 freddo, 1caldo, file altri)\n - nome init file\n - seed (<0)\n - numero misure\n - passi decoerenza\n - passi termalizzazione");
        return 0;
    }
    //valori default parametri
    FILE *init_f = NULL, *output;
    long seed_v = -42, *seed = &seed_v;
    double beta = 0.3;
    int Nterm = 0,      //numero di passi da saltare per termalizzare (moltiplicato per L^2)
        Ndecoer = 1,    //numero di passi da saltare tra due misure (moltiplicato per L^2)
        Nmis = 100000,  //numero di misure
        init_flag = 1,  //1 = caldo, 0 = freddo
        Nlato = 16;     //lato del reticolo
    //lettura parametri
    switch(argc) {
        case 10: sscanf(argv[9], "%d", &Nterm);
        case 9: sscanf(argv[8], "%d", &Ndecoer);
        case 8: sscanf(argv[7], "%d", &Nmis);
        case 7: sscanf(argv[6], "%d", seed);
        case 6:
        case 5: if (argv[4][0] == '0' && argv[4][1] == '\0') init_flag = 0;
                else if (argv[4][0] == '1' && argv[4][1] == '\0') init_flag = 1;
                else init_f = fopen(argv[4], "r");
        case 4: sscanf(argv[3], "%d", &Nlato);
        case 3: sscanf(argv[2], "%lf", &beta);
        case 2: output = fopen(argv[1], "a"); break;
        case 1: output = fopen("Default_out.txt", "w"); break;
        default: perror("Troppi argomenti");
            output = fopen("Default_out.txt", "w");
    }
    //printf("%d\n", Nmis);
    //init
    int *field = malloc(sizeof(int)*Nlato*Nlato);
    if (init_f == NULL) {
        if (init_flag) init_caldo(field, seed, Nlato);
        else init_freddo(field, Nlato);
    }
    else {
        init_file(field, Nlato, init_f, seed);
        fclose(init_f);
    }
    int geo[Nlato][Ngeo]; geometry(geo, Nlato);
    fprintf(output, "#%f\t%d\t%d\t%d\t%d\n", beta, Nlato, Nmis, Nterm, Ndecoer);
    //campionamento
    Ndecoer *= Nlato*Nlato;
    Nterm *= Nlato*Nlato;
    for (int i=0; i<Nterm; i++) {
        passo_metropolis(field, Nlato, beta, seed, geo);
    }
    for (int i=0; i<Nmis; i++) {
        for (int j=0; j<Ndecoer; j++){
            passo_metropolis(field, Nlato, beta, seed, geo);
        }
        fprintf(output, "%d\t%d\n", energia(field, Nlato, geo), tot_magnetizza(field, Nlato));
    }
    //salvataggio config
    if (argc == 5) init_f = fopen(argv[4], "w");
    else if (argc >= 6) init_f = fopen(argv[5], "w");
    else init_f = fopen("LastConfig.txt", "w");
    fprintf(init_f, "%d\n", seed_v);
    print_field(field, Nlato, init_f);
    free(field);
    fclose(output); fclose(init_f);
    return 0;
}
