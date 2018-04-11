import MySQLdb

usuario = ""
senha = ""
base = ""
host = ""

db = MySQLdb.connect(host, usuario, senha, base)
cursor = db.cursor()

sql_registro = "SELECT * FROM registro ORDER BY id DESC"

cursor.execute(sql_registro)
res_registro = cursor.fetchall()

print("ID       n_proc      n_passos        temp(incr)        rede        contorno        duracao     data")

for row in res_registro:
    id = row[0]
    n_proc = row[1]
    n_passos = row[2]
    temp_i = row[3]
    temp_f = row[4]
    temp_incr = row[5]
    nx = row[6]
    ny = row[7]
    nz = row[8]
    cx = row[9]
    cy = row[10]
    duracao = row[11]
    data_inicio = row[12]
    cz = row[13]
    
    print("%d       %d      %d      %f-%f(%f)       %dx%dx%d        %d%d%d      %f      %s" %(id, n_proc, n_passos, temp_i, \
    temp_f, temp_incr, nx, ny, nz, cx, cy, cz, duracao, data_inicio))

id_dados = int(input("\n\nDigite a id desejada: "))

sql_registro_filtro = "SELECT * FROM registro WHERE id=%d" %id_dados 
cursor.execute(sql_registro_filtro)
res_registro_filtro = cursor.fetchall()

for row in res_registro_filtro:
    id = row[0]
    n_proc = row[1]
    n_passos = row[2]
    temp_i = row[3]
    temp_f = row[4]
    temp_incr = row[5]
    nx = row[6]
    ny = row[7]
    nz = row[8]
    cx = row[9]
    cy = row[10]
    duracao = row[11]
    data_inicio = row[12]
    cz = row[13] 

nome_arquivo = "dados_%d_%d_[%dx%dx%d].dat" %(id, n_passos, nx, ny, nz) 
arquivo = open(nome_arquivo, "w")

sql_dados = "SELECT * FROM dados where id_registro=%d ORDER BY temp ASC" %(id_dados)
cursor.execute(sql_dados)
res_dados = cursor.fetchall()

nome_arquivo = "dados_%d_%d_[%dx%dx%d].dat" %(id, n_passos, nx, ny, nz) 
arquivo = open(nome_arquivo, "w")

arquivo.write("ID: %d\n" %(id))
arquivo.write("No processadores: %d\n" %(n_proc))
arquivo.write("No passos: %d\n" %(n_passos))
arquivo.write("Temperatura: [%f, %f], incremento: %f\n" %(temp_i, temp_f, temp_incr))
arquivo.write("Rede: %dx%dx%d     Condicao de contorno: %dx%dx%d\n" %(nx, ny, nz, cx, cy, cz))
arquivo.write("Duracao: %f    Data: %s\n\n" %(duracao, data_inicio))

for row in res_dados:  
    arquivo.write("%.3f %.20f %.20f %.20f %.20f %.20f %.20f %.20f\n" %(row[2], row[3], row[4], row[5], row[6], row[7], row[8], row[9]))

print("\n-- Arquivo salvo! --\n")
arquivo.close()
db.close()