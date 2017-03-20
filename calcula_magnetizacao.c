/*

	Descrição: Esta função calcula a magnetizacao em uma
		       simulação do modelo de Ising em 2d.
	Autor: Rodrigo Carvalho (rpatriciocarvalho@gmail.com)
	Última modificação: 16/02/2016

	Observações:
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
