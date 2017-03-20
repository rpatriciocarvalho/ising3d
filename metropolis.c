/*

	Descrição: Esta função executa o algoritmo de Metropolis em uma
		       simulação do modelo de Ising em 2d.
	Autor: Rodrigo Carvalho (rpatriciocarvalho@gmail.com)
	Última modificação: 6/11/2014

	Observações:
*/

#include "funcoes_ising2d.h"
#include <math.h>
#include <stdlib.h>

void metropolis(float temperatura){

    int x, y, z;
    double exponencial;
    int novo_spin, soma_vizinhos, delta_energia;

    for(x = 0; x < NX; x++){
        for(y = 0; y < NY; y++){
            for(z = 0; z < NZ; z++){
                if(VIZINHO == 0) {
                    soma_vizinhos = vizinho_nulo(x,y,z,1) +
                                    vizinho_nulo(x,y,z,2) +
                                    vizinho_nulo(x,y,z,3) +
                                    vizinho_nulo(x,y,z,4) +
                                    vizinho_nulo(x,y,z,5) +
                                    vizinho_nulo(x,y,z,6);
                } else {
                    soma_vizinhos = vizinho_periodico(x,y,z,1) +
                                    vizinho_periodico(x,y,z,2) +
                                    vizinho_periodico(x,y,z,3) +
                                    vizinho_periodico(x,y,z,4) +
                                    vizinho_periodico(x,y,z,5) +
                                    vizinho_periodico(x,y,z,6);
                }

                novo_spin = -1*rede[x][y][z];
                delta_energia = -J*(novo_spin - rede[x][y][z])*soma_vizinhos;

                if (delta_energia < 0){
                    rede[x][y][z] *= -1;
                } else {
                    exponencial = (double) pow(M_E, (double) (-delta_energia)/(K_B*temperatura));
                    if(exponencial > ranmar()) rede[x][y][z] *= -1;
                }
            }
        }
    }
}
