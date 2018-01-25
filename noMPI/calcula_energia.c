/*

	Descrição: Esta função calcula a energia em uma
		   simulação do modelo de Ising.
	Autor: Rodrigo Carvalho (rpatriciocarvalho@gmail.com)
*/

#include "funcoes.h"

double calcula_energia(){

    int x, y, z;
    double energia = 0;

    /*
    A função vizinho possui 4 parâmetros:
            1o - posição x do sítio analisado;
            2o - posição y do sítio analisado;
            3o - vizinho (1 = direito, 2 = superior, 3 = esquerdo, 4 = inferior)
    */

    for(x = 0; x < NX; x++){
        for(y = 0; y < NY; y++){
            for(z = 0; z < NZ; z++){
                
                    energia += rede[x][y][z]*(vizinho(x,y,z,1, VIZINHO_X) +
                        vizinho(x,y,z,2, VIZINHO_Y) +
                        vizinho(x,y,z,3, VIZINHO_X) +
                        vizinho(x,y,z,4, VIZINHO_Y) +
                        vizinho(x,y,z,5, VIZINHO_Z) +
                        vizinho(x,y,z,6, VIZINHO_Z));

            }
        }
    }

    energia *= -J;
    return(energia);
}
