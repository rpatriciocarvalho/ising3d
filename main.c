/*

	Descri��o: Simula��o do modelo de Ising em 3d.
	Autor: Rodrigo Carvalho (rpatriciocarvalho@gmail.com)
	�ltima modifica��o: 23/01/2016

	Observa��es: Os par�metros da simula��o devem se modificados
				no arquivo 'funcoes_ising2d.h'
	14			
				25/01/2017 - Adicionei o calculo do cumulante de Binder
*/

#include <stdio.h>
#include <math.h>
#include <time.h>
#include <stdlib.h>
#include "funcoes_ising2d.h"

int main(){

	int i, passos_termalizacao, quantidade_medidas;
	double magnetizacao, energia, temperatura, fator_normalizacao;
	double calor_especifico, suscetibilidade_mag, mag, ene;
    double magnetizacao_2, energia_2, magnetizacao_3, magnetizacao_4;
	double desvio_magnetizacao, desvio_quarto_magnetizacao, desvio_energia;
    double cumulante;
	char nome_arquivo[25];

	// Iniciando o gerador de n�meros aleat�rios
    srand(time(NULL)); // Inicia semente
    rmarin((int) rand()%31328, (int) rand()%30081); //inicia o ranmar

	/*
	Se "TERMALIZACAO" for diferente de zero, ent�o
	gera o gr�fico da energia por passos e n�o executa
	a simula��o.
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
	Executa a simula��o caso n�o tenhamos que
	calcular a termaliza��o.
	*/

	if (TERMALIZACAO == 0) {
        for(quantidade_medidas = 1; quantidade_medidas <= MEDIDAS; quantidade_medidas++){
            FILE *arquivo_dados;

            if (CLUSTER == 0){
                if(VIZINHO == 0) {
                    sprintf(nome_arquivo,"dados_%dx%dx%d_%d_[%d]_nulo.dat", NX, NY, NZ, N_PASSOS, quantidade_medidas);
                } else {
                    sprintf(nome_arquivo,"dados_%dx%dx%d_%d_[%d]_periodico.dat", NX, NY, NZ, N_PASSOS, quantidade_medidas);
                }
                arquivo_dados = fopen(nome_arquivo, "w");
            }

            passos_termalizacao = N_PASSOS*0.1; // Definindo 10% do n�mero total de passos para termalizar.
            fator_normalizacao = (double) NX*NY*NZ*(N_PASSOS-passos_termalizacao+1.0);

            for(temperatura=TEMP_I;temperatura<=TEMP_F; temperatura+=INCRE_TEMP){

                magnetizacao = energia = 0.0;
                magnetizacao_2 = magnetizacao_3 = magnetizacao_4 = energia_2 = 0.0;
                
                for(i=1; i<=N_PASSOS; i++){

                    metropolis(temperatura);

                     if (i >= passos_termalizacao){
                        mag = fabsf(calcula_magnetizacao());
                        ene = fabsf(calcula_energia());
                        magnetizacao += mag;
                        energia += ene;
                        magnetizacao_2 += pow(mag, 2);
                        magnetizacao_3 += pow(mag, 3);
                        magnetizacao_4 += pow(mag, 4);
                        energia_2 += pow(ene, 2);
                    }
                }

                magnetizacao /= fator_normalizacao;
                energia /= fator_normalizacao;
                magnetizacao_2 /= fator_normalizacao*NX*NY*NZ;
                magnetizacao_3 /= fator_normalizacao*pow(NX*NY*NZ, 2);
                magnetizacao_4 /= fator_normalizacao*pow(NX*NY*NZ, 3);
                energia_2 /= fator_normalizacao*NX*NY*NZ;

                calor_especifico = (double) (NX*NY*NZ)*(energia_2 - pow(energia, 2))/(K_B*pow(temperatura, 2));
                suscetibilidade_mag = (double) (NX*NY*NZ)*(magnetizacao_2 - pow(magnetizacao, 2))/(K_B*temperatura);
                desvio_magnetizacao = (double) sqrt(magnetizacao_2 - pow(magnetizacao, 2));
                desvio_quarto_magnetizacao = (double) magnetizacao_4 - 4*magnetizacao_3*magnetizacao +
                                                6*magnetizacao_2*pow(magnetizacao, 2) - 3*pow(magnetizacao, 4);
                desvio_energia = (double) sqrt(energia_2 - pow(energia, 2));
                cumulante = 1-(desvio_quarto_magnetizacao/(3*pow(magnetizacao_2 - pow(magnetizacao, 2), 2)));

				if(temperatura >= 3.5 && temperatura <= 5.0){
					cumulante = cumulante;
				} else {
					cumulante = 0;
				}
					
				
                // Ser� escrito no arquivo as seguintes colunas
                // Temp. -  Mag. m�dia - desvio mag.
                // - energia m�dia - desvio energia - calor espec?fico - suscetibilidade_mag - cumulante

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
        fclose(arquivo_dados);
        }
	}
	return 0;
}
