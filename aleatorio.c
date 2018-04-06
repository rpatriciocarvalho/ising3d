/*

	Descrição: Esta função gera um número aleatório com o RANMAR em uma
           simulação do modelo de Ising.
    Autor: Rodrigo Carvalho (rpatriciocarvalho@gmail.com)
*/

#include "funcoes.h"
#include "mt19937ar.h"

int aleatorio(){

    // Esta função gera um número aleatório 1 e -1

    if(genrand_real2() < 0.5){
        return(1);
    } else {
        return(-1);
    }

}