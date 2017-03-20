/*

	Descri��o: Esta fun��o calcula a energia em uma
		       simula��o do modelo de Ising em 2d.
	Autor: Rodrigo Carvalho (rpatriciocarvalho@gmail.com)
	�ltima modifica��o: 16/02/2016

	Observa��es:
*/

#include "funcoes_ising2d.h"

double calcula_energia(){

	int x, y, z;
	double energia = 0;

	/*
	A fun��o vizinho possui 4 par�metros:
		1o - posi��o x do s�tio analisado;
		2o - posi��o y do s�tio analisado;
		3o - vizinho (1 = direito, 2 = superior, 3 = esquerdo, 4 = inferior)
	*/

    for(x = 0; x < NX; x++){
		for(y = 0; y < NY; y++){
			for(z = 0; z < NZ; z++){
                if(VIZINHO == 0) {
                    energia += rede[x][y][z]*(vizinho_nulo(x,y,z,1) +
									vizinho_nulo(x,y,z,2) +
									vizinho_nulo(x,y,z,3) +
									vizinho_nulo(x,y,z,4) +
									vizinho_nulo(x,y,z,5) +
									vizinho_nulo(x,y,z,6));
                } else {
                    energia += rede[x][y][z]*(vizinho_periodico(x,y,z,1) +
									vizinho_periodico(x,y,z,2) +
									vizinho_periodico(x,y,z,3) +
									vizinho_periodico(x,y,z,4) +
									vizinho_periodico(x,y,z,5) +
									vizinho_periodico(x,y,z,6));
                }
			}
		}
	}

	energia *= -J;
	return(energia);
}
