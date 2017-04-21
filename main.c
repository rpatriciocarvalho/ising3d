/*
    Descrição: Simulação do modelo de Ising em 3d.
    Autor: Rodrigo Carvalho (rpatriciocarvalho@gmail.com)
    Última modificação: 20/03/2017
    Observações: Os parâmetros da simulação devem se modificados
                 no arquivo 'funcoes_ising2d.h'
                 25/01/2017 - Adicionei o calculo do cumulante de Binder
*/

#include <stdio.h>
#include <math.h>
#include <time.h>
#include <stdlib.h>
#include <mpi.h>
#include "funcoes_ising2d.h"

int main(int argc, char** argv){

    int i, passos_termalizacao, quantidade_medidas;
    int n_passos_mpi, inicio_node, fim_node, fator_normalizacao;
    double magnetizacao, energia, temperatura;
    double calor_especifico, suscetibilidade_mag, mag, ene;
    double magnetizacao_2, energia_2, magnetizacao_4;
    double desvio_magnetizacao, desvio_energia;
    double cumulante;
    double magnetizacao_node, energia_node, magnetizacao_2_node, magnetizacao_4_node, energia_2_node;
    char nome_arquivo[25];

    time_t inicio, fim;

    inicio = time(NULL);
    
    // Iniciando o gerador de números aleatórios
    srand(time(NULL)); // Inicia semente
    rmarin((int) rand()%31328, (int) rand()%30081); //inicia o ranmar

    /*
    Se "TERMALIZACAO" for diferente de zero, então
    gera o gráfico da energia por passos e não executa
    a simulação.
    */

    iniciar_rede(); // Gerando a rede inicial

    if (CLUSTER == 0){
        if (TERMALIZACAO != 0){

            fator_normalizacao = (double) NX*NY;
            FILE *fp = fopen("termalizacao.dat", "w");

            for(i=1; i <= N_PASSOS; i++){
                metropolis(1.8); // 1.8
                energia = calcula_energia()/fator_normalizacao;
                magnetizacao = calcula_magnetizacao()/fator_normalizacao;
                fprintf(fp, "%d %f %f\n", i, energia, magnetizacao);
            }
            fclose(fp);
        }
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
            
            n_passos_mpi = N_PASSOS/nprocs;

            for(i=0; i<nprocs; i++){
                if(myrank == i){
                        inicio_node = i*n_passos_mpi + 1;
                        fim_node = (i+1)*n_passos_mpi;
                        passos_termalizacao = i*n_passos_mpi + n_passos_mpi*0.1;                    
                }
            }

            fator_normalizacao = (int) (N_PASSOS - N_PASSOS*0.1)*NX*NY*NZ;
                
                FILE *arquivo_dados;
                
                if (CLUSTER == 0){
                    if(VIZINHO == 0) {
                        sprintf(nome_arquivo,"dados_%dx%dx%d_%d_[%d]_nulo.dat", NX, NY, NZ, N_PASSOS, quantidade_medidas);
                    } else if (VIZINHO == 1) {
                        sprintf(nome_arquivo,"dados_%dx%dx%d_%d_[%d]_unitario.dat", NX, NY, NZ, N_PASSOS, quantidade_medidas);
                    } else {
                        sprintf(nome_arquivo,"dados_%dx%dx%d_%d_[%d]_periodico.dat", NX, NY, NZ, N_PASSOS, quantidade_medidas);
                    }
                    arquivo_dados = fopen(nome_arquivo, "w");
                }
            
            for(temperatura=TEMP_I;temperatura<=TEMP_F; temperatura+=INCRE_TEMP){
                
                iniciar_rede();
                magnetizacao = energia = 0.0;
                magnetizacao_2 = magnetizacao_4 = energia_2 = 0.0;
                magnetizacao_node = energia_node = magnetizacao_2_node = 0.0;
                magnetizacao_4_node = energia_2_node = 0.0;
                
                for(i=inicio_node; i <= fim_node; i++){
                    
                    metropolis(temperatura);

                     
                     if (i > passos_termalizacao){
                        mag = fabsf(calcula_magnetizacao());
                        ene = fabsf(calcula_energia());
                        magnetizacao_node += mag;
                        energia_node += ene;
                        magnetizacao_2_node += pow(mag, 2);
                        magnetizacao_4_node += pow(mag, 4);
                        energia_2_node += pow(ene, 2);
                        
                        MPI_Reduce(&magnetizacao_node, &magnetizacao, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
                        MPI_Reduce(&energia_node, &energia, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
                        MPI_Reduce(&magnetizacao_2_node, &magnetizacao_2, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
                        MPI_Reduce(&magnetizacao_4_node, &magnetizacao_4, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
                        MPI_Reduce(&energia_2_node, &energia_2, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
                       
                    }
                      
                }

                if (myrank == 0) {
                    
                    magnetizacao /= (double) fator_normalizacao;
                    energia /= (double) fator_normalizacao;
                    magnetizacao_2 /= (double) fator_normalizacao*NX*NY*NZ;
                    magnetizacao_4 /= (double) fator_normalizacao*pow(NX*NY*NZ, 3);
                    energia_2 /= (double) fator_normalizacao*NX*NY*NZ;

                    calor_especifico = (double) (NX*NY*NZ)*(energia_2 - pow(energia, 2))/(K_B*pow(temperatura, 2));
                    suscetibilidade_mag = (double) (NX*NY*NZ)*(magnetizacao_2 - pow(magnetizacao, 2))/(K_B*temperatura);
                    desvio_magnetizacao = (double) sqrt(magnetizacao_2 - pow(magnetizacao, 2));
                    desvio_energia = (double) sqrt(energia_2 - pow(energia, 2));
                    cumulante = 1-(magnetizacao_4/(3*pow(magnetizacao_2, 2)));
                    
                    // Será escrito no arquivo as seguintes colunas
                    // Temp. -  Mag. média - desvio mag.
                    // - energia média - desvio energia - calor específico - suscetibilidade_mag - cumulante

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
                        printf("%f %2.20f %2.20f %2.20f %2.20f %2.20f %2.20f\n", temperatura,
                                    magnetizacao,
                                    desvio_magnetizacao,
                                    energia,
                                    desvio_energia,
                                    calor_especifico,
                                    suscetibilidade_mag);
                    }
                }
                
            }
        fclose(arquivo_dados);
        MPI_Finalize();
        }
    }

    fim = time(NULL);
    printf("O tempo de execucao em segundos é %f\n", difftime(fim, inicio));
    return 0;
}