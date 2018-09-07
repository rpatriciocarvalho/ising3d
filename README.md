# ising3d
Simulação do modelo de Ising 3D programado em C para computação paralela usando openMPI e base de dados MySQL.

Trabalhamos com o método de Monte Carlo, algoritmo de Metropolis e o gerador de números aletórios ranmar.

Há dois arquivos importantes na pasta `ferramentas`. O primeiro, `criaTabelas.sql` refere-se a criação das tabelas na base de dados MySQL (ou MariaDB) que irão armazenar os dados produzidos pela simulação. 

O segundo, `recuperaDB.py`, é uma simples aplicação em python que verifica quais foram as simulação realizadas e em seguida, ao escolher um dos registros, salva em um arquivo `.dat` os dados gerados pela simulação.

Grandezas calculadas na simulação:
* Temperatura
* Magnetização
* Desvio-padrão da magnetização
* Energia
* Desvio-padrão da energia
* Susceptibilidade Magnética
* Calor específico
* Cumulante de Binder
