/*

	Descri��o: Esta fun��o gera a rede inicial em uma
		       simula��o do modelo de Ising em 2d.
	Autor: Rodrigo Carvalho (rpatriciocarvalho@gmail.com)
	�ltima modifica��o: 16/02/2016

	Observa��es:
*/

#include "funcoes_ising2d.h"

void iniciar_rede(){

    int x, y, z;

	/*
	Gera rede inicial. Se "PARTIDA" for nulo, ent�o �
	iniciada como todos os sitios unit�rios. Se for
	n�o nulo, ent�o � aleat�ria a rede inicial.
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
