/*

	Descrição: Esta função calcula a magnetizacao em uma
		       simulação do modelo de Ising.
	Autor: Rodrigo Carvalho (rpatriciocarvalho@gmail.com)

*/

#include "funcoes_ising2d.h"

double calcula_magnetizacao(){

    int x, y, z, magnetizacao;

    magnetizacao = 0;

    for(x = 0; x < NX; x++){
        for(y = 0; y < NY; y++){
            for(z = 0; z < NZ; z++){
                magnetizacao += rede[x][y][z];
            }
        }
    }

    return (magnetizacao);

}
