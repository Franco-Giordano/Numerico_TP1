# K.u = F


def crearK(n):
	"""Crea la matriz K de coeficientes"""

	if n < 5:
		return [[0.0]*5 for x in range(5)]

	dim = n+1
	K = [[0.0]*dim for x in range(dim)]

	K[0][0] = 1.0
	K[-1][-1] = 1.0

	vectorFijo = [-4.0, 5.0, -4.0, 1.0]

	for i, p in enumerate(vectorFijo):
		K[1][i] = p
		K[-2][-i-1] = p

	vectorFijo2 = [1.0] + vectorFijo
	vectorFijo2[2] += 1

	for i in range(2,n-1):
		offset = i-2
		for j,e in enumerate(vectorFijo2):
			K[i][j+offset] = e

	return K

def q(x, g=6.0):
	return g + g**2 * (x - x**2) 

def x(n, i, L=1.0):
	return ((i*L) / n)

def f(n,i,k=1.0,L=1.0):
	xi = x(n,i)

	return (q(xi) * (L/n)**4) / k

def crearF(n):
	"""Crear matriz F"""

	dim = n+1

	F = [0.0]*dim

	for i in range(1,n):
		F[i] = f(n,i)

	return F


def SOR(K, F, anterior, w=1.0):

	actualGS = GaussSeidel(K, F, anterior)

	actual = [anterior[i] * (1.0 - w) + w * actualGS[i] for i in range(len(anterior))]

	return actual


def GaussSeidel(K, F, anterior):

	result = anterior[::]

	for i in range(len(K)):

		factor = 0.0
		for j in range(len(K)):
			if j != i:
				factor += (K[i][j]*result[j]) / K[i][i]

		result[i] = (F[i]/K[i][i]) - factor

	return result

def normaInfinito(a):
	return abs(max(a, key=lambda x: abs(x)))

def criterioConvergencia(actual, anterior, tolerancia=0.01):
	resta = [actual[i]-anterior[i] for i in range(len(actual))]

	normaResta = normaInfinito(resta)

	normaActual = normaInfinito(actual)

	return normaResta/normaActual <= tolerancia

def main():

	F = crearF(5)
	K = crearK(5)
	ant = [0.0,1.0,0.0,10.0,0.0,0.0]

	actual = SOR(K,F,ant, w=1.03)

	while not (criterioConvergencia(actual, ant)):
		ant = actual
		actual = SOR(K,F,ant,w=1.03)
		print(actual)


K = [[10.0,2.0,6.0], [1.0,10.0,4.0], [2.0,-7.0,-10.0]]
F = [28.0,7.0,-17.0]
ant = [1.0,2.0,3.0]

for _ in  range(8):
	ant = SOR(K,F,ant, w=1.033)
	print(ant)

#main()