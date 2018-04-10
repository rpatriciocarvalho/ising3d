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
#include <mpi.h>
#include <unistd.h>
#include "mt19937ar.h"
#include "funcoes.h"
#include <mysql/mysql.h>
#include <string.h>

int main(int argc, char** argv){

    unsigned long int i, passos_termalizacao, k;
    double magnetizacao, energia, temperatura, fator_normalizacao;
    double calor_especifico, suscetibilidade_mag, mag, ene;
    double magnetizacao_2, energia_2, magnetizacao_4;
    double desvio_magnetizacao, desvio_energia;
    double cumulante;

    // Variáveis para computação paralela
    //double m_temp_mpi;
    int inicio_node, fim_node;

    time_t inicio; // Variáveis para marcar o tempo de execução
	time_t fim;

    inicio = time(NULL); // Início da execução da aplicação
        
    // Iniciando o gerador de números aleatórios
    srand(time(NULL)); // Inicia semente
    init_genrand((int) rand()%31328);

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
        int id_registro = 0;
        char query_dados[2200];
            
        // Conexao com a base de dados
        MYSQL conexao;
        mysql_init(&conexao);
        mysql_real_connect(&conexao, HOST, DBUSUARIO, DBSENHA, BASEDADOS, 0, NULL, 0);

        // Inclusão do MPI
        int myrank, nprocs;

        MPI_Init(&argc, &argv);
        MPI_Comm_size(MPI_COMM_WORLD, &nprocs);
        MPI_Comm_rank(MPI_COMM_WORLD, &myrank);
        
        if (myrank == 0){
            // Inserindo dados de registro
            char query_registro[2000];
            
            sprintf(query_registro, "INSERT INTO registro (n_proc, n_passos, temp_i, temp_f, temp_incr, nx, ny, nz, cx, cy, cz, duracao) "
                "values(%d, %d, %f, %f, %f, %d, %d, %d, %d, %d, %d, 0)",  
                nprocs, N_PASSOS, TEMP_I, TEMP_F, INCRE_TEMP, NX, NY, NZ, VIZINHO_X, VIZINHO_Y, VIZINHO_Z);
            mysql_query(&conexao, query_registro);
            
            // Capturando ID da última query
            id_registro = mysql_insert_id(&conexao);
        }

        // Enviando a id da última quere para todas as unidades de processamento
        MPI_Bcast(&id_registro, 1, MPI_INT, 0, MPI_COMM_WORLD);

        int tamanho_array_temperaturas_medidas = 0;
        tamanho_array_temperaturas_medidas = (TEMP_F-TEMP_I)/INCRE_TEMP;

        float temperaturas_medidas[tamanho_array_temperaturas_medidas];

        for(i = 0; i < tamanho_array_temperaturas_medidas; i++){
            temperaturas_medidas[i] = TEMP_I + i*INCRE_TEMP;
        }

        inicio_node = myrank*(tamanho_array_temperaturas_medidas/nprocs);
        fim_node = inicio_node + (tamanho_array_temperaturas_medidas/nprocs) - 1;
        
        /*
        // Divide a temperatura por unidade de processamento
        m_temp_mpi = (double) (TEMP_F*1.0 - TEMP_I*1.0)/(nprocs*1.0);

        // Estabelece temperaturas iniciais e finais para cada unidade de processamento
        inicio_node = TEMP_I + (myrank*1.0)*m_temp_mpi + INCRE_TEMP*(1.0);
        fim_node = inicio_node + m_temp_mpi; 
        */

        // Definindo 10% do número total de passos para termalizar.
        passos_termalizacao = N_PASSOS*0.1;
        
        // Fator de normalização que será usado nas médias
        fator_normalizacao = (unsigned long int) (N_PASSOS - N_PASSOS*0.1)*NX*NY*NZ;

        // É calculado as medidas para cada temperatura
        for(k=inicio_node; k <= fim_node; k++){
            
            // "Limpando" a variável query_dados
            memset(&query_dados, 0, sizeof(query_dados));
            
            // Cada unidade de processamento iniciará após um determinado tempo
            sleep((int) myrank + 1); 
            
            // A temperatura que sera calculada
            temperatura = temperaturas_medidas[k];

            iniciar_rede(); 
            
            // Zerando variáveis
            magnetizacao = energia = 0.0;
            magnetizacao_2 = magnetizacao_4 = energia_2 = 0.0;

            // Variáveis referentes a porcentagem do programa
            double j = 1.0;
            double porcentagem_passos;

            // Para cada unidade de processamento os passos de Monte Carlo
            for(i=1; i <= N_PASSOS; i++){
                
                metropolis(temperatura); // Algoritmo de Metropolis
                
                // Começa a medir apenas depois da termalização
                if (i > passos_termalizacao){
                    mag = fabsf(calcula_magnetizacao());
                    ene = fabsf(calcula_energia());
                    magnetizacao += mag;
                    energia += ene;
                    // Variáveis para cálculo dos desvios e outras grandezas
                    magnetizacao_2 += pow(mag, 2);
                    magnetizacao_4 += pow(mag, 4);
                    energia_2 += pow(ene, 2);                        
                }

                // Porcentagem do progresso da simulação                  
                porcentagem_passos = ((i*1.0)/N_PASSOS)*100.0;

                if(porcentagem_passos >= j) {
                    if (myrank == 0){   
                        printf("%2.2f%% - %2.2f%%\n", ( (temperatura-inicio_node)/(fim_node-inicio_node)*100, porcentagem_passos);
                        j++;
                    }   
                }                    
            }

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

            // Inserindo dados de registro
            sprintf(query_dados, "INSERT INTO dados (id_registro, temp, mag, desv_mag, energia, desv_energia, calor_esp, sus_mag, cumulante) "
                "values(%d, %f, %2.20f, %2.20f, %2.20f, %2.20f, %2.20f, %2.20f, %2.20f)",  
                id_registro, temperatura, magnetizacao, desvio_magnetizacao, energia, desvio_energia, calor_especifico, suscetibilidade_mag, cumulante);
            mysql_query(&conexao, query_dados);
            
            /*              
                Escrevemos os dados em um arquivo.
                
                Será escrito no arquivo as seguintes colunas:
                    Temperatura -  Magnetização média - desvio magnetização.
                    energia média - desvio energia - calor específico 
                    suscetibilidade magnética - cumulante de Binder
            
           
            fprintf(arquivo_dados, "%f %2.20f %2.20f %2.20f %2.20f %2.20f %2.20f %2.20f\n", temperatura,
                                        magnetizacao,
                                        desvio_magnetizacao,
                                        energia,
                                        desvio_energia,
                                        calor_especifico,
                                        suscetibilidade_mag,
                                        cumulante);

            */
        }
        
        //fclose(arquivo_dados);

        if (myrank == 0){
            // Registra o momento que o processo termina e imprime na tela
            fim = time(NULL);
            printf("Tempo: %f\n", difftime(fim, inicio));
            
            // Inserindo dados de registro
            char query_registro_update[100];
            sprintf(query_registro_update, "UPDATE registro SET duracao=%f WHERE id=%d", difftime(fim, inicio), id_registro);
            mysql_query(&conexao, query_registro_update);
        }

        MPI_Finalize(); // Finalizamos o MPI        

        mysql_close(&conexao);    
    }

    return 0;
}
