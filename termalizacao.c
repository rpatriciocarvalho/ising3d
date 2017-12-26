#include "funcoes_ising2d.h"
#include <stdio.h>
#include <math.h>

void calcula_termalizacao(){
    
    int i;
    double magnetizacao, energia, temperatura, fator_normalizacao;
    
    fator_normalizacao = (double) NX*NY*NY;
    FILE *fp = fopen("termalizacao.dat", "w");

    for(i=1; i <= N_PASSOS; i++){
        metropolis(1.8); // Temperatura 1.8
        energia = fabsf(calcula_energia()/fator_normalizacao);
        magnetizacao = fabsf(calcula_magnetizacao()/fator_normalizacao);
        fprintf(fp, "%d %f %f\n", i, energia, magnetizacao);
    }
    fclose(fp);

}