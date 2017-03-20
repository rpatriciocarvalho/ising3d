/*

	Descri��o: Este cabe�alho possui os prot�tipos de fun��es da
		       simula��o do modelo de Ising em 2d.
	Autor: Rodrigo Carvalho (rpatriciocarvalho@gmail.com)
	�ltima modifica��o: 23/01/2017

	Observa��es:
*/

// Par�metros ---------------------------
#define N_PASSOS	100000	// N�mero de passos de Monte Carlo
#define NX		8	// Dimens�o x da rede
#define NY		8	// Dimens�o y da rede
#define NZ              8	// Dimens�o z da rede
#define TEMP_I		0.1	// Temperatura inicial da rede
#define TEMP_F		10	// Temperatura final da rede
#define INCRE_TEMP	0.1	// Incremento da temperatura
#define J		1	// Constante de acoplamento
#define K_B		1	// Constate de Boltzmann
#define	PARTIDA		0	// Partida fria (=0) ou quente (!=0)
#define TERMALIZACAO	0	// Verifica a termaliza��o (!=0) ou n�o (=0)
#define CLUSTER		0	// 1 => Simula��o no cluster; 0 => Simula��o no pc
#define MEDIDAS     1   // N�mero de medidas que ser�o feitas
#define VIZINHO     1 // 0 = vizinho nulo; 1 = vizinho periodico

double calcula_energia();
double calcula_magnetizacao();
void metropolis(float temperatura);
void iniciar_rede();
int aleatorio();
int vizinho_nulo(int x, int y, int z, int posicao);
int vizinho_periodico(int x, int y, int z, int posicao);

// Criando a rede
int rede[NX][NY][NZ];

// Gerador de numeros aleat�rios Ranmar --------------
double ranmar(void);
void rmarin(int, int);

static double u[97], c, cd, cm;
int i97, j97;
//----------------------------------------------------
