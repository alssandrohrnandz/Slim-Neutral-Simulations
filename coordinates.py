# -*- coding: utf-8 -*-

import pandas as pd
import os
import math
import sys

# Verificar que se ha pasado un argumento
if len(sys.argv) != 2:
    print("Uso: python coordenadas.py <SGE_TASK_ID>")
    sys.exit(1)

# Obtener el SGE_TASK_ID del argumento
SGE_TASK_ID = int(sys.argv[1])

# Tamaño de la matriz
matrix_size = 100

# Verificar que matrix_size es un cuadrado perfecto
n = int(math.sqrt(matrix_size))
if n * n != matrix_size:
    raise ValueError("matrix_size debe ser un cuadrado perfecto")

# Función para obtener Long y Lat dinámicamente
def get_coordinates(pob, n):
    # Extraer el índice de la población (p1 -> 1, p2 -> 2, ..., pn -> n)
    index = int(pob[1:]) - 1
    # Calcular Long y Lat basados en la matriz nxn
    row = index // n
    col = index % n
    Long = n - row
    Lat = col + 1
    return Long, Lat

carpeta = "ArchivosListos"

# Cargar el archivo txt en un DataFrame
archivo = f'sim_{SGE_TASK_ID}.txt'
if not os.path.exists(archivo):
    print(f"El archivo {archivo} no existe.")
    sys.exit(1)

df = pd.read_csv(archivo, delimiter=' ', header=None, names=['Date_mean','AlleleFrequency', 'Pob', 'AlleleCount','DerivedAlleles'])

# Aplicar las funciones a las columnas 'Pob' y 'Long', y almacenar los resultados en nuevas columnas
df[['Long', 'Lat']] = df['Pob'].apply(lambda x: pd.Series(get_coordinates(x, n)))

# Guardar el DataFrame procesado en un nuevo archivo txt
output_file = os.path.join(carpeta, f'SimPro_{SGE_TASK_ID}.txt')
df.to_csv(output_file, sep='\t', index=False)
print(f"Archivo procesado y guardado en {output_file}")
