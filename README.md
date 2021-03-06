# FlujoCancerMama
## Participantes:
* Patricia Trujillo Rodríguez (Patriciatr)
* Irene Sánchez Jiménez (IreneSanx)
* Lucía Valverde Martínez (LuciaVM)
* Clara Jiménez Valverde (clarajvv)

Sobre el flujo:
En este flujo de datos podremos realizar diferentes DEAs (Análisis de Expresión Diferencial), en los que siempre enfrentaremos a la variable de supervivencia del paciente con la que elija el usuario (en algunos casos también podrá elegir los rangos o valores dentro de la misma) de las disponibles. Estas serán:
1. Subtipos de tumor, entre los que encontramos Luminal A, Basal y TNBC (triple negativo). Podrá elegir cualquier par de ellas.
2. Número de nodos afectados, que pueden ser ninguno, uno, dos, o tres. El usuario podrá elegir cualquier par de ellos.
3. Edad del paciente. En este caso el usuario podrá elegir los rangos a comparar.
4. Grado del tumor. El grado del tumor podrá tener 3 grados posibles. Pero se comprobará la supervivencia entre el grado T1 y T3. Siendo estos: T1 (de grado bajo), las células parecen ser como las células ordinarias y son de desarrollo moderado; T3 (de alto grado), las células parecen ser extremadamente diferentes a las células típicas y se están desarrollando más rápidamente.
5. Si existen nodos afectados. A diferencia de la segunda opción, medimos solamente si hay algún nodo afectado o no. Es el equivalente de comparar los cero nodos con el resto de variables. Al ser los valores posibles positivo o negativo, el usuario no tiene la opción de elegirlos.
6. Estado ER, referido a si el tumor posee receptores de estrógenos. Los valores posibles son de nuevo positivo o negativo, por lo que el usuario no los elegirá.
7. Estado PR, referido a si el tumor posee receptores de progesterona. Como en el caso anterior, los valores serán positivo o negativo, por lo que el usuario no los elegirá.

Para la ejecución del script de bash, debe hacerse desde la carpeta FlujoCancerMama que se descargará al clonar el repositorio. Al llamar al script sin ningún parámetro, se arrojará un mensaje con la información sobre qué parametros hay que pasarle en función del análisis que se desee hacer. 
