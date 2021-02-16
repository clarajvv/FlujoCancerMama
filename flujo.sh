#/bin/bash
if [ "$1" == "" ];then 
	echo "En este flujo de datos podremos realizar diferentes DEAs, enfrentando siempre la variable de supervivencia del paciente con la que elija el usuario. 

Puede elegir entre las siguientes opciones (cuyo número pasará como primer parámetro):

1. Suptipos de tumor. Tenemos Basal (1), LuminalA (2) y TNBC (3). Tendrá que elegir un par de ellas, pasando los números asociados como segundo y tercer parámetro. 

2. Número de nodos afectados, que pueden ser ninguno (0), uno (1), dos (2) o tres (3). Debe elegir un par de ellos, pasando sus números asociados como segundo y tercer parámetro.

3. Edad del paciente. Debe elegir como segundo y tercer parámetro las edades para separar el conjunto. La edad mínima es 26 y la máxima 90.

4. Grado del tumor. Aquí se comprobará la supervivencia entre el grado T1 y T3, por lo que no tiene que introducir ningún parámetro más. 

5. Si existen o no nodos afectados. Tampoco tiene que introducir ningún parámetro adicional.

6. Estado ER, referido a si el tumor posee receptores de estrógenos. No necesita ningún parámetro más.

7. Estado PR, sobre si el tumor tiene receptores de progesterona. No requiere ningún parámetro extra."

else

Rscript Code/Settings/1_getData.R
if [ "$1" == "1" ]; then
	Rscript Code/Options/2_1_SubsettingTumourSubtypes.R $2 $3
	Rscript Code/Options/3_DEA.R
elif [ "$1" == "2" ]; then
	Rscript Code/Options/2_2_SubsettingNodeNumber.R $2 $3
	Rscript Code/Options/3_DEA.R
elif [ "$1" == "3" ]; then
	Rscript Code/Options/2_3_subsettingAges.R $2 $3
	Rscript Code/Options/3_DEA.R
elif [ "$1" == "4" ]; then
	Rscript Code/Options/2_4_subsettingGrade.R
	Rscript Code/Options/3_DEA.R
elif [ "$1" == "5" ]; then
	Rscript Code/Options/2_5_SubsettingNodes.R
	Rscript Code/Options/3_DEA.R
elif [ "$1" == "6" ]; then
	Rscript Code/Options/2_6_SubsettingERStatus.R
	Rscript Code/Options/3_DEA.R
elif [ "$1" == "7" ]; then
	Rscript Code/Options/2_7_SubsettingPRStatus.R
	Rscript Code/Options/3_DEA.R
else
	echo "El número introducido no corresponde con ninguna de las opciones disponibles"
fi

fi