#include "funcoes.h"

int vizinho(int x, int y, int z, int posicao, int condicao){
    
    if(condicao == 0){
        vizinho_nulo(x, y, z, posicao);
    } else {
        vizinho_periodico(x, y, z, posicao);
    }

}