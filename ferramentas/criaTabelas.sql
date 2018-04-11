CREATE TABLE dados (
`id` INT NOT NULL AUTO_INCREMENT PRIMARY KEY,
id_registro int not null,
temp float not null,
mag double not null,
desv_mag double not null,
energia double not null,
desv_energia double not null,
calor_esp double not null,
sus_mag double not null,
cumulante double not null
);

CREATE TABLE registro (
`id` INT NOT NULL AUTO_INCREMENT PRIMARY KEY,
n_proc int not null,
n_passos int not null,
temp_i float not null,
temp_f float not null,
temp_incr float not null,
nx int not null,
ny int not null,
nz int not null,
cx int not null,
cy int not null,
duracao float not null,
data_inicio TIMESTAMP DEFAULT CURRENT_TIMESTAMP ON UPDATE CURRENT_TIMESTAMP,
cz int not null
);