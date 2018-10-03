def crearMatrizK(n):

	"""Crea la matriz K de coeficientes provenientes de la ecuacion diferencial discretizada."""

	if n < 5:
		return [[0.0]*5 for x in range(5)]

	dimension = n+1

	K = [[0.0] * dimension for x in range(dimension)]

	K[0][0] = 1.0
	K[-1][-1] = 1.0

	#Como notamos un patron en las filas de la matriz K
	#Adopatmos un vector fijo que se va trasladando y modificando

	vectorFijo = [-4.0, 5.0, -4.0, 1.0]

	for i, p in enumerate(vectorFijo):

		K[1][i] = p
		K[-2][-i-1] = p

	vectorFijo2 = [1.0] + vectorFijo
	vectorFijo2[2] += 1

	for i in range(2, n-1):

		offset = i - 2

		for j, e in enumerate(vectorFijo2):

			K[i][j + offset] = e

	return K

def q(x, grupo = 6.0):

	"""Calcula la carga, donde el parametro "grupo" es estatico."""

	#"grupo" se puede pasar opcionalmente.

	return grupo + grupo ** 2 * (x - x ** 2) 

def x(n, i, L = 1.0):

	"""Calcula el argumento de la funcion "q"."""

	return ((i * L) / n)

def f(n, i, k = 1.0, L = 1.0):

	xi = x(n,i)

	return (q(xi) * (L / n) ** 4) / k

def crearF(n):

	"""Crear matriz F"""

	dimension = n + 1

	F = [0.0] * dimension

	for i in range(1,n):

		F[i] = f(n,i)

	return F


def SOR(K, F, anterior, w = 1.0):

	actualGS = GaussSeidel(K, F, anterior)

	actual = [anterior[i] * (1.0 - w) + w * actualGS[i] for i in range(len(anterior))]

	return actual

def newSOR(K, F, anterior, w =1.0):

	res = anterior[::]

	for j in range(len(res)):
		posGS = singleGaussSeidel(K, F, res, j)

		res[j] = (posGS - res[j])*w + res[j]

	return res

def GaussSeidel(K, F, anterior):

	resultado = anterior[::]

	for i in range(len(K)):

		factor = 0.0
		for j in range(len(K)):
			if j != i:
				factor += (K[i][j] * resultado[j]) / K[i][i]

		resultado[i] = (F[i] / K[i][i]) - factor

	return resultado

def singleGaussSeidel(K, F, anterior, pos):

	resultado = anterior[::]

	factor = 0.0

	for j in range(len(K)):
		if j != pos:
			factor += (K[pos][j] * resultado[j]) / K[pos][pos]

		resultado[pos] = (F[pos] / K[pos][pos]) - factor

	return resultado[pos]

def normaInfinito(vector):

	return abs(max(vector, key = lambda x: abs(x)))

def criterioConvergencia(actual, anterior, tolerancia = 0.01):

	resta = [actual[i] - anterior[i] for i in range(len(actual))]

	normaResta = normaInfinito(resta)

	normaActual = normaInfinito(actual)

	return normaResta / normaActual <= tolerancia

def dumpLista(lista, dumpFile):

	"""dump de los datos de la lista seleccionada, formateado para mejor exportacion"""

	for posicionElemento in range(len(lista)):

		if posicionElemento == 0:
			dumpFile.write("%.04f" % lista[posicionElemento])

		else:
			dumpFile.write("	" + "%.04f" % lista[posicionElemento])

	dumpFile.write("\n")

def dumpDatosGrafico(iteracionesTotalesPorW, factoresDeRelajacion, dumpFile):

	"""Dump de los datos necesarios para graficar"""

	for posicionElemento in range(len(iteracionesTotalesPorW)):
		dumpFile.write("%d    %.2f\n" % (iteracionesTotalesPorW[posicionElemento], factoresDeRelajacion[posicionElemento]))


def main():

	tomarIntervalos = True
	intervalos = []

	while tomarIntervalos:

		aux = int(input("Ingrese un tamanio de intervalo mayor a 4, caso contrario finaliza lectura: "))

		if aux > 4:
			intervalos.append(aux)

		else:
			tomarIntervalos = False

	if intervalos:
		print("Se ejecutara el programa para los siguientes intervalos: {} ".format(intervalos))

	factoresDeRelajacion = [x / 100 for x in range(100, 200, 5)]
	iteracionesTotalesPorW = []

	for n in intervalos:

		dumpFile = open("Intervalos_" + str(n) + ".txt", "w")

		for factorDeRelajacion in factoresDeRelajacion:

			dimension = n + 1

			anterior = [0] * dimension

			K = crearMatrizK(n)

			F = crearF(n)

			actual = newSOR(K, F, anterior, factorDeRelajacion)

			iteracionesTotales = 1

			dumpFile.write("Procesamiento del SEL con factor de relajacion = " + str(factorDeRelajacion) + "\n\n")
			dumpLista(actual, dumpFile)

			while not (criterioConvergencia(actual, anterior)):

				anterior = actual
				actual = newSOR(K, F, anterior, factorDeRelajacion)

				dumpLista(actual, dumpFile)

				iteracionesTotales += 1

			iteracionesTotalesPorW.append(iteracionesTotales)

			dumpFile.write("\n" + "Iteraciones totales: " + str(iteracionesTotales) + "\n\n")
			dumpFile.write("----------------------------------------------------- \n\n")

		dumpDatosGrafico(iteracionesTotalesPorW, factoresDeRelajacion, dumpFile)

		dumpFile.close()


"""K = [[10.0,2.0,6.0], [1.0,10.0,4.0], [2.0,-7.0,-10.0]]
F = [28.0,7.0,-17.0]
ant = [1.0,2.0,3.0]

for _ in  range(8):
	ant = SOR(K,F,ant, w=1.033)
	print(ant)"""

main()