/*

	Descrição: Esta função encontra o vizinho em uma
		       simulação do modelo de Ising.
	Autor: Rodrigo Carvalho (rpatriciocarvalho@gmail.com)

*/

#include "funcoes.h"

int vizinho_periodico(int x, int y, int z, int posicao){

    int eixo_x, eixo_y, eixo_z;

    if (posicao == 1){
        eixo_x = x + 1;
        eixo_y = y;
        eixo_z = z;
        if (eixo_x == NX) eixo_x = 0;
    }

    if (posicao == 2){
        eixo_x = x;
        eixo_y = y - 1;
        eixo_z = z;
        if (eixo_y == -1) eixo_y = (NY - 1);
    }

    if (posicao == 3){
        eixo_x = x - 1;
        eixo_y = y;
        eixo_z = z;
        if (eixo_x == -1) eixo_x = (NX - 1);
    }

    if (posicao == 4){
        eixo_x = x;
        eixo_y = y + 1;
        eixo_z = z;
        if (eixo_y == NY) eixo_y = 0;
    }

    if (posicao == 5){
        eixo_x = x;
        eixo_y = y;
        eixo_z = z - 1;
        if (eixo_z == -1) eixo_z = (NZ - 1);
    }

    if (posicao == 6){
        eixo_x = x;
        eixo_y = y;
        eixo_z = z + 1;
        if (eixo_z == NZ) eixo_z = 0;
    }

    return (rede[eixo_x][eixo_y][eixo_z]);
}