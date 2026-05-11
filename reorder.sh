#!/bin/bash

# 1. Definir el directorio base
DIRECTORIO_BASE="Out/linear"

# 2. Inicializar un array vacío
FOLDER_LIST=()

# 3. Llenar el array SOLO con subcarpetas
# El /*/ al final asegura que solo busque directorios
for path in "$DIRECTORIO_BASE"/*/; do
    
    # Esta comprobación evita errores si el directorio base está vacío
    [ -d "$path" ] || continue 
    
    # Usamos basename para quitar la ruta (Out/circle/) y quedarnos solo con el nombre de la carpeta
    NOMBRE_CARPETA=$(basename "$path")
    IFS='_' read -r -a CHAR_ARRAY <<< "$NOMBRE_CARPETA"
    N=${CHAR_ARRAY[0]:2}
    E=${CHAR_ARRAY[1]:2}
    Temp=${CHAR_ARRAY[2]:5}
    mu=${CHAR_ARRAY[3]:3}
    gamma=${CHAR_ARRAY[4]:2}
    N_pot=$(echo "l($N)/l(2)" | bc -l)
    # We do the same process to find out all the subfolders
    for subpath in "$DIRECTORIO_BASE/$NOMBRE_CARPETA"/*/; do
        SUB_NOMBRE_CARPETA=$(basename "$subpath")
        IFS='_' read -r -a SUBCHAR_ARRAY <<< "$SUB_NOMBRE_CARPETA"
        R=${SUBCHAR_ARRAY[0]:6}
        nT=${SUBCHAR_ARRAY[1]:3}
        measT=${SUBCHAR_ARRAY[2]:6}
        stT=${SUBCHAR_ARRAY[3]:4}
        M=$(echo "scale=0; sqrt($N)" | bc)
        #echo "$R $nT $measT $stT $M"
        NEW_FILE_NAME="G=${gamma}_E=${E}_Temp=${Temp}_mu=${mu}/N=${N_pot:0:2}_M=${M}_R=${R}_nT=${nT}_measT=${measT}_stT=${stT}"
        echo "$DIRECTORIO_BASE/$NOMBRE_CARPETA/$SUB_NOMBRE_CARPETA"
        echo "$DIRECTORIO_BASE/$NEW_FILE_NAME"
        # Copy last file into file with new data org
        mkdir -p "$DIRECTORIO_BASE/$NEW_FILE_NAME"
        cp -r "$DIRECTORIO_BASE/$NOMBRE_CARPETA/$SUB_NOMBRE_CARPETA/." "$DIRECTORIO_BASE/$NEW_FILE_NAME"
    done
done
