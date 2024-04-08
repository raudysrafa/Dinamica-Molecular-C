# Dinamica-Molecular-C
En esta aplicación se simula un gas de Argón interactuando a través de un potencial de Lennard-Jones, 
se calculan las fuerzas que actúan sobre cada partícula, asi como la aceleración y velocidad en cada instante,
se guardan en un fichero los resultados de los cálculos de las energías cinética y potencial.
Se hace un uso extensivo del trabajo con arreglos, se implementa una tabla de Verlet que es una rutina que rastrea la posición
de las partículas más cercanas a la partícula examinada. Esto evita realizar un bucle sobre todas las partículas de la simulación
con la respectiva ganancia de tiempo
