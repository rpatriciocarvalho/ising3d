#include "funcoes.h"
#include <stdio.h>
#include <math.h>

void calcula_termalizacao(){
    
    int i, j;
    double magnetizacao, energia, temperatura, fator_normalizacao;
    double porcentagem_passos;

    fator_normalizacao = (double) NX*NY*NY;
    FILE *fp = fopen("termalizacao.dat", "w");

    j = 1.0;

    for(i=1; i <= N_PASSOS; i++){
        metropolis(4.5);
        energia = fabsf(calcula_energia()/fator_normalizacao);
        magnetizacao = fabsf(calcula_magnetizacao()/fator_normalizacao);
        fprintf(fp, "%d %f %f\n", i, energia, magnetizacao);

        porcentagem_passos = ((i*1.0)/N_PASSOS)*100.0;

        if(porcentagem_passos >= j) {
            printf("%2.0f%%\n", porcentagem_passos);
            j++;   
        }
    }
    fclose(fp);

}