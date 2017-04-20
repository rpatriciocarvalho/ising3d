﻿/*

	Descriçãoo: Este cabeçalho possui os protótipos de funções da
		       simulação do modelo de Ising em 2d.
	Autor: Rodrigo Carvalho (rpatriciocarvalho@gmail.com)
	Última modificação: 23/01/2017

	Observações:
*/

// Parâmetros ---------------------------
#define N_PASSOS            10000	// Número de passos de Monte Carlo
#define NX                  8           // Dimensão x da rede
#define NY                  8           // Dimensão y da rede
#define NZ                  8           // Dimensão z da rede
#define TEMP_I              1         // Temperatura inicial da rede
#define TEMP_F              8          // Temperatura final da rede
#define INCRE_TEMP          0.1         // Incremento da temperatura
#define J                   1           // Constante de acoplamento
#define K_B                 1           // Constate de Boltzmann
#define	PARTIDA             0           // Partida fria (=0) ou quente (!=0)
#define TERMALIZACAO        0           // Verifica a termalização (!=0) ou não (=0)
#define CLUSTER             0           // 1 => Simulação no cluster; 0 => Simulação no pc
#define MEDIDAS             1           // Número de medidas que serão feitas
#define VIZINHO             2           // 0 = vizinho nulo; 1 = vizinho unitario; > 1 = periodico

double calcula_energia();
double calcula_magnetizacao();
void metropolis(float temperatura);
void iniciar_rede();
int aleatorio();
int vizinho_nulo(int x, int y, int z, int posicao);
int vizinho_periodico(int x, int y, int z, int posicao);
int vizinho_unitario(int x, int y, int z, int posicao);

// Criando a rede
int rede[NX][NY][NZ];

// Gerador de numeros aleatórios Ranmar --------------
double ranmar(void);
void rmarin(int, int);

static double u[97], c, cd, cm;
int i97, j97;