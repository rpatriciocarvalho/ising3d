/*

	Descrição: Esta função gera a rede inicial em uma
		       simulação do modelo de Ising.
	Autor: Rodrigo Carvalho (rpatriciocarvalho@gmail.com)

*/

#include "funcoes.h"

void iniciar_rede(){

    int x, y, z;

    /*
    Gera rede inicial. Se "PARTIDA" for nulo, então é
    iniciada como todos os sitios unitários. Se for
    não nulo, então é aleatória a rede inicial.
    */

    if (PARTIDA == 0) {
        for(x=0; x < NX; x++){
            for(y=0; y < NY; y++){
                for(z = 0; z < NZ; z++){
                        rede[x][y][z] = 1;
                }
            }
        }
    } else {
        for(x=0; x < NX; x++){
            for(y=0; y < NY; y++){
                for(z = 0; z < NZ; z++){
                        rede[x][y][z] = aleatorio();
                }
            }
        }
    }
}