/*

    Descrição: Simulação do modelo de Ising em 3d.
    Autor: Rodrigo Carvalho (rpatriciocarvalho@gmail.com)

    Observações: Os parâmetros da simulação devem se modificados
                 no arquivo 'funcoes_ising2d.h'

*/

#include <stdio.h>
#include <math.h>
#include <time.h>
#include <stdlib.h>
#include <mpi.h>
#include "funcoes_ising2d.h"

int main(int argc, char** argv){

    unsigned long int i, passos_termalizacao, quantidade_medidas;
    double magnetizacao, energia, temperatura, fator_normalizacao;
    double calor_especifico, suscetibilidade_mag, mag, ene;
    double magnetizacao_2, energia_2, magnetizacao_4;
    double desvio_magnetizacao, desvio_energia;
    double cumulante;
    char nome_arquivo[25];

    // Variáveis para computação paralela
    double m_temp_mpi;
    double inicio_node, fim_node, temperatura_node;
    double magnetizacao_node, energia_node, magnetizacao_2_node;
    double magnetizacao_4_node, energia_2_node;
    double calor_especifico_node, suscetibilidade_mag_node, desvio_magnetizacao_node;
    double desvio_energia_node, cumulante_node;


    time_t inicio; // Variáveis para marcar o tempo de execução
	time_t fim;

    inicio = time(NULL); // Início da execução da aplicação
        
    // Iniciando o gerador de números aleatórios
    srand(time(NULL)); // Inicia semente
    rmarin((int) rand()%31328, (int) rand()%30081); //inicia o ranmar

   /*
    Se "TERMALIZACAO" for diferente de zero, então
    gera o gráfico da energia por passos e não executa
    a simulação.

    Apenas calcula a termalização.
    */
    
    if (TERMALIZACAO != 0){
        iniciar_rede(); // Gerando a rede inicial
        calcula_termalizacao();
    } 

    /*
    Executa a simulação caso não tenhamos que
    calcular a termalização.
    */

    if (TERMALIZACAO == 0) {
        for(quantidade_medidas = 1; quantidade_medidas <= MEDIDAS; quantidade_medidas++){
            
            // Inclusão do MPI
            int myrank, nprocs;

            MPI_Init(&argc, &argv);
            MPI_Comm_size(MPI_COMM_WORLD, &nprocs);
            MPI_Comm_rank(MPI_COMM_WORLD, &myrank);
            
            m_temp_mpi = (double) TEMP_F/(nprocs - 1); // Divide a temperatura por unidade de processamento

            // Estabelece temperaturas iniciais e finais para cada unidade de processamento
            for(i=1; i<nprocs; i++){
                if(myrank == i){
                    if (myrank == 1){ 
                        inicio_node = INCRE_TEMP;
                        fim_node = m_temp_mpi;
                    } else {
                        inicio_node = (i-1)*m_temp_mpi;
                        fim_node = inicio_node + m_temp_mpi;
                    }                    
                }
            }

            // Definindo 10% do número total de passos para termalizar.
            passos_termalizacao = N_PASSOS*0.1;
            
            // Fator de normalização que será usado nas médias
            fator_normalizacao = (unsigned long int) (N_PASSOS - N_PASSOS*0.1)*NX*NY*NZ;

            // Cria um arquivo para armazenar os dados    
            FILE *arquivo_dados;
                        
            sprintf(nome_arquivo,"dados_%dx%dx%d_%u_[%d]_[%d-%d-%d].dat", NX, NY, NZ, N_PASSOS, 
                quantidade_medidas, VIZINHO_X, VIZINHO_Y, VIZINHO_Z);

            if(myrank == 0) {
                arquivo_dados = fopen(nome_arquivo, "w"); // O arquivo é criado apenas no master
                temperatura = 0;
            }
            
            // Zerando variáveis
            magnetizacao = energia = 0.0;
            magnetizacao_2 = magnetizacao_4 = energia_2 = 0.0;

            if(myrank != 0) {
                // É calculado as medidas para cada temperatura
                for(temperatura_node=inicio_node;temperatura_node<=fim_node; temperatura_node+=INCRE_TEMP){
                
                    iniciar_rede(); 
                    
                    // Zerando variáveis
                    magnetizacao_node = energia_node = magnetizacao_2_node = 0.0;
                    magnetizacao_4_node = energia_2_node = 0.0;
                    
                    // Para cada unidade de processamento os passos de Monte Carlo
                    for(i=1; i <= N_PASSOS; i++){
                        
                        metropolis(temperatura_node); // Algoritmo de Metropolis
                        
                        // Começa a medir apenas depois da termalização
                        if (i > passos_termalizacao){
                            mag = fabsf(calcula_magnetizacao());
                            ene = fabsf(calcula_energia());
                            magnetizacao_node += mag;
                            energia_node += ene;
                            // Variáveis para cálculo dos desvios e outras grandezas
                            magnetizacao_2_node += pow(mag, 2);
                            magnetizacao_4_node += pow(mag, 4);
                            energia_2_node += pow(ene, 2);                        
                        }                    
                    }

                    // Envia os dados de cada unidade de processamento para a unidade de processamento primária
                    MPI_Gather(&temperatura_node, 1, MPI_DOUBLE, &temperatura, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
                    MPI_Gather(&magnetizacao_node, 1, MPI_DOUBLE, &magnetizacao, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
                    MPI_Gather(&energia_node, 1, MPI_DOUBLE, &energia, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
                    MPI_Gather(&magnetizacao_2_node, 1, MPI_DOUBLE, &magnetizacao_2, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
                    MPI_Gather(&magnetizacao_4_node, 1, MPI_DOUBLE, &magnetizacao_4, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
                    MPI_Gather(&energia_2_node, 1, MPI_DOUBLE, &energia_2, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
                }
            }
            
            // Se o processo estiver rodando na unidade de processamento primária, calcula-se as médias
            if (myrank == 0) {
                    
                magnetizacao /= (double) fator_normalizacao;
                energia /= (double) fator_normalizacao;
                magnetizacao_2 /= (double) fator_normalizacao*NX*NY*NZ;
                magnetizacao_4 /= (double) fator_normalizacao*pow(NX*NY*NZ, 3);
                energia_2 /= (double) fator_normalizacao*NX*NY*NZ;
                
                // Calculamos outras grandezas
                calor_especifico = (double) (NX*NY*NZ)*(energia_2 - pow(energia, 2))/(K_B*pow(temperatura, 2));
                suscetibilidade_mag = (double) (NX*NY*NZ)*(magnetizacao_2 - pow(magnetizacao, 2))/(K_B*temperatura);
                desvio_magnetizacao = (double) sqrt(magnetizacao_2 - pow(magnetizacao, 2));
                desvio_energia = (double) sqrt(energia_2 - pow(energia, 2));
                cumulante = 1-(magnetizacao_4/(3*pow(magnetizacao_2, 2)));
                
                /*              
                    Escrevemos os dados em um arquivo.
                    
                    Será escrito no arquivo as seguintes colunas:
                        Temperatura -  Magnetização média - desvio magnetização.
                        energia média - desvio energia - calor específico 
                        suscetibilidade magnética - cumulante de Binder
                */

                
                printf("%2.2f%%\n", (float) temperatura); // Porcentagem da simulacao
                fprintf(arquivo_dados, "%f %2.20f %2.20f %2.20f %2.20f %2.20f %2.20f %2.20f\n", temperatura,
                                    magnetizacao,
                                    desvio_magnetizacao,
                                    energia,
                                    desvio_energia,
                                    calor_especifico,
                                    suscetibilidade_mag,
                                    cumulante);                
            }
        
            if(myrank == 0) fclose(arquivo_dados); // Fechamos o arquivo
            MPI_Finalize(); // Finalizamos o MPI        
        }
    }

    // Registra o momento que o processo termina e imprime na tela
    fim = time(NULL);
    printf("O tempo de execucao em segundos foi %f\n", difftime(fim, inicio));
    
    return 0;
}