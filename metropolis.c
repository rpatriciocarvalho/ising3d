/*

	Descrição: Esta função executa o algoritmo de Metropolis em uma
		       simulação do modelo de Ising em 3d.
	Autor: Rodrigo Carvalho (rpatriciocarvalho@gmail.com)

*/

#include "funcoes.h"
#include <math.h>
#include <stdlib.h>

void metropolis(float temperatura){

    int x, y, z;
    double exponencial;
    int novo_spin, soma_vizinhos, delta_energia;

    for(x = 0; x < NX; x++){
        for(y = 0; y < NY; y++){
            for(z = 0; z < NZ; z++){
                
                soma_vizinhos = vizinho(x,y,z,1, VIZINHO_X) +
                                   vizinho(x,y,z,2, VIZINHO_Y) +
                                   vizinho(x,y,z,3, VIZINHO_X) +
                                   vizinho(x,y,z,4, VIZINHO_Y) +
                                   vizinho(x,y,z,5, VIZINHO_Z) +
                                   vizinho(x,y,z,6, VIZINHO_Z);
                
                novo_spin = -1*rede[x][y][z];
                delta_energia = -(J/(K_B*temperatura))*(novo_spin - rede[x][y][z])*soma_vizinhos;

                if (delta_energia < 0 || log(ranmar()) < -delta_energia ){
                    rede[x][y][z] *= -1;
                }
            }
        }
    }
}