/*

    Descrição: Simulação do modelo de Ising em 3d.
    Autor: Rodrigo Carvalho (rpatriciocarvalho@gmail.com)

    Observações: Os parâmetros da simulação devem se modificados
                 no arquivo 'funcoes.h'

*/

#include <stdio.h>
#include <math.h>
#include <time.h>
#include <stdlib.h>
#include "funcoes.h"

int main(int argc, char** argv){

    int i, passos_termalizacao, quantidade_medidas;
    double magnetizacao, energia, temperatura, fator_normalizacao;
    double calor_especifico, suscetibilidade_mag, mag, ene;
    double magnetizacao_2, energia_2, magnetizacao_4;
    double desvio_magnetizacao, desvio_energia;
    double cumulante;
    char nome_arquivo[25];

    time_t inicio, fim; // Variáveis para marcar o tempo de execução

    inicio = time(NULL); // Início da execução da aplicação
    
    // Iniciando o gerador de números aleatórios
    srand(time(NULL)); // Inicia semente
    rmarin((int) rand()%31328, (int) rand()%30081); //inicia o ranmar

    /*
    Se "TERMALIZACAO" for diferente de zero, então
    gera o gráfico da energia por passos e não executa
    a simulação.

    Apenas calcula a termalização se não for em um cluster.
    */
    
    if (TERMALIZACAO != 0 && CLUSTER == 0){
        iniciar_rede(); // Gerando a rede inicial
        calcula_termalizacao();
    } 

    /*
    Executa a simulação caso não tenhamos que
    calcular a termalização.
    */

    if (TERMALIZACAO == 0) {
        // Podemos fazer várias rodadas da aplicação em uma única execução
        for(quantidade_medidas = 1; quantidade_medidas <= MEDIDAS; quantidade_medidas++){
            
            // Cria-se o arquivo de dados
            FILE *arquivo_dados; 
                
            if (CLUSTER == 0){
                    
                sprintf(nome_arquivo,"dados_%dx%dx%d_%d_[%d]_[%d-%d-%d].dat", NX, NY, NZ, N_PASSOS, 
                        quantidade_medidas, VIZINHO_X, VIZINHO_Y, VIZINHO_Z);
                arquivo_dados = fopen(nome_arquivo, "w");
            }
            
            // Para cada temperatura dentro do intervalo
            for(temperatura=TEMP_I;temperatura<=TEMP_F; temperatura+=INCRE_TEMP){
                
                // Iniciamos a rede
                iniciar_rede();
            
                // Zeramos as variáveis
                magnetizacao = energia = 0.0;
                magnetizacao_2 = magnetizacao_4 = energia_2 = 0.0;                  
                
                // Definindo 10% do número total de passos para termalizar.
                passos_termalizacao = N_PASSOS*0.1;

                // Estabelecemos o fator de normalização e média
                fator_normalizacao = (double) NX*NY*NZ*(N_PASSOS-passos_termalizacao+1.0);
              
                // Os passos do método de Monte Carlo
                for(i=0; i <= N_PASSOS; i++){
                    
                    // Execução do algoritmo de Metropolis para dada temperatura
                    metropolis(temperatura);

                    /* 
                    Apenas contabilizamos para a média os os passos
                    acima da termalização
                    */

                    if (i >= passos_termalizacao){
                        mag = fabsf(calcula_magnetizacao()); // O módulo da magnetização
                        ene = fabsf(calcula_energia()); // O módulo da energia
                        
                        /*
                        Somamos as magnetizações e energias para cada passo,
                        assim como o quadrado e o quarta potência para cálculos futuros
                        */

                        magnetizacao += mag; 
                        energia += ene;
                        magnetizacao_2 += pow(mag, 2);
                        magnetizacao_4 += pow(mag, 4);
                        energia_2 += pow(ene, 2);     
                    }
                }

                // Normalizamos as medidas e calculamos as médias
                magnetizacao /= fator_normalizacao;
                energia /= fator_normalizacao;
                magnetizacao_2 /= fator_normalizacao*NX*NY*NZ;
                magnetizacao_4 /= fator_normalizacao*pow(NX*NY*NZ, 3);
                energia_2 /= fator_normalizacao*NX*NY*NZ;

                // Calculamos outras grandezas
                calor_especifico = (double) (NX*NY*NZ)*(energia_2 - pow(energia, 2))/(K_B*pow(temperatura, 2));
                suscetibilidade_mag = (double) (NX*NY*NZ)*(magnetizacao_2 - pow(magnetizacao, 2))/(K_B*temperatura);
                desvio_magnetizacao = (double) sqrt(magnetizacao_2 - pow(magnetizacao, 2));
                desvio_energia = (double) sqrt(energia_2 - pow(energia, 2));
                cumulante = 1-(magnetizacao_4/(3*pow(magnetizacao_2, 2)));
                    
                /*              
                Se a simulação não estiver rodando em um cluster, escrevemos os dados em um arquivo.
                Caso contrário escrevemos na tela.

                Será escrito no arquivo as seguintes colunas:
                    Temperatura -  Magnetização média - desvio magnetização.
                    energia média - desvio energia - calor específico 
                    suscetibilidade magnética - cumulante de Binder
                */

                if (CLUSTER == 0){
                    printf("[%d/%d] - %2.2f%%\n", quantidade_medidas, MEDIDAS, (float) (temperatura/TEMP_F)*100.0); // Porcentagem da simulacao
                    fprintf(arquivo_dados, "%f %2.20f %2.20f %2.20f %2.20f %2.20f %2.20f %2.20f\n", temperatura,
                                            magnetizacao,
                                            desvio_magnetizacao,
                                            energia,
                                            desvio_energia,
                                            calor_especifico,
                                            suscetibilidade_mag,
                                            cumulante);
                } else {
                    printf("%f %2.20f %2.20f %2.20f %2.20f %2.20f %2.20f %2.20f\n", temperatura,
                                magnetizacao,
                                desvio_magnetizacao,
                                energia,
                                desvio_energia,
                                calor_especifico,
                                suscetibilidade_mag,
                                cumulante);
                }                                
            }
        
            // Fechamos o arquivo
            fclose(arquivo_dados);
 
        }
    }

    // Imprimimos na tela o tempo de execução da aplicação
    fim = time(NULL);
    printf("O tempo de execucao em segundos %f\n", difftime(fim, inicio));
    return 0;
}