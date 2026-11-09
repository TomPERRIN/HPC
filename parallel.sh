
#!/bin/bash

# Nombre de points par simulation
NB_POINTS=1000000000

# Nombre de réplications
N=10

# Fichier de résultats
OUTPUT="resultats.txt"

# Initialisation du fichier de résultats
> "$OUTPUT"

# Lancement des simulations en parallèle
for i in $(seq 1 $N)
do
    SEED=$((1000 + i))

    echo "Lancement de l'experience $i avec seed $SEED"

    ./simPi "$NB_POINTS" "$SEED" >> "$OUTPUT" &

done

# Attente de la fin de tous les processus
wait

echo "Toutes les simulations sont terminees."

echo "Resultats :"
cat "$OUTPUT"
