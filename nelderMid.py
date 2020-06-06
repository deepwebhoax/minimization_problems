import matplotlib.pyplot as plt
from numpy import array
import random
import numpy as np
from math import sqrt

#f = lambda x: 4*(x[0]-3)**2 + (x[1]-2)**2
# Функция Розенброка
a, b = 1, 100
f = lambda x: (1-x[0])**2 + 100*(x[1]-x[0]**2)**2
# sliding decent mod (while)
# h = [lambda x: 25-x[0]**2-x[1]**2]
g=[]
h = [lambda x: x[0]**2 + x[1]**2 - 9* x[1] + 4.25]
def nelderMid(phi, f=lambda x: (a-x[0])**2 + b*(x[1]-x[0]**2)**2, xStart=[5,5], simpSizeCoef = 0.01):
	n = 3
	kvcf = 0
	alpha, beta, gamma = 1, (-0.5, 0.5), 2
	# Начальные точки
	xStart = array(xStart)
	x_f = ( xStart-array([0, simpSizeCoef]), None),( xStart+array([simpSizeCoef, 0]), None),( xStart+array([0, simpSizeCoef]) , None)
	
	x_f = [(p[0], f(p[0])) for p in x_f]
	kvcf += 3
	project = lambda xh, xc, teta: xh + (1+teta)*(xc-xh)

	previousF = ['*','*','*', '*', '*']
	k=0
	i=0
	thing = 0
	# превычисление (для первой итерации)
	x_f.sort(key=lambda p:p[1])
	xNew = x_f[0][0]
	while phi - f(x_f[0][0]) < 0:
		x_f.sort(key=lambda p:p[1])
		xc = (x_f[0][0]+x_f[1][0])/2

		teta = 1
		xNew = project(x_f[2][0], xc, teta)
		fNew = f(xNew)
		kvcf +=1
		if fNew<x_f[0][1]:
			teta = gamma
		elif fNew<x_f[1][1]:
			teta = alpha
		elif fNew<x_f[2][1]:
			teta = beta[1]
		else:
			teta = beta[0]
			x1 = project(x_f[0][0], x_f[1][0], teta)
			x2 = project(x_f[0][0], x_f[2][0], teta)
			x_f = [(x1,f(x1)), (x2, f(x2)), x_f[0]]
			x_f.sort(key=lambda p:p[1])
			i += 1
			kvcf += 2
			continue
			
		xNew = project(x_f[2][0], xc, teta)
		fNew = f(xNew)
		kvcf += 1
		# Проверка на зацикливание
		if fNew == previousF[-2] or fNew == previousF[-3] or fNew == previousF[-4]:
			# print("in")
			teta=beta[0]
			xNew = project(x_f[2][0], xc, teta)
			fNew = f(xNew)
			kvcf += 1
			
		x_f[2] = (xNew, fNew)
		i+=1
		if i >5000:
			return None
		previousF.append(fNew)
		x_f.sort(key=lambda p:p[1])
	x_f.sort(key=lambda p:p[1])	
	return x_f[0][0], kvcf


# Функція основного методу
def nelderMid1(f=lambda x: (a-x[0])**2 + b*(x[1]-x[0]**2)**2, xStart=[5,5], startSize=0.5, e = 10e-6):
	n = 3
	e = 10e-6
	kvcf = 0
	alpha, beta, gamma = 1, (-0.5, 0.5), 2
	# Начальные точки
	x_f = (xStart-array([startSize/2, 0]), None),(xStart+array([startSize/2, 0]), None),(xStart+array([0, sqrt(startSize**2-(startSize/2)**2)]) , None)
	# x_f = (xStart-array([0, 0.1**simpSizeCoef]), None),(xStart+array([0.1**simpSizeCoef, 0]), None),(xStart+array([0, 0.1**simpSizeCoef]) , None)
	l = [x_f[0][0], x_f[1][0], x_f[2][0]]
	kvcf +=3
	x_f = [(p[0], f(p[0])) for p in x_f]
	xc = (x_f[0][0]+x_f[1][0])/2
	# print("Начальные точки: "+" ".join(["x"+str(i+1)+" = "+str(x_f[i][0]) for i in range(3)]))

	# print(" ".join(["f"+str(i+1)+" = "+str(x_f[i][1]) for i in range(n)]))

	project = lambda xh, xc, teta: xh + (1+teta)*(xc-xh)

	previousF = ['*','*','*', '*', '*']
	k=0
	i=0
	# превычисление (для первой итерации)
	x_f.sort(key=lambda p:p[1])
	while (sum([(x_f[i][1]-f(xc))**2 for i in range(n)])/(n+1))**(0.5) > e:
		# print(str(i+1)+")")
		x_f.sort(key=lambda p:p[1])
		xc = (x_f[0][0]+x_f[1][0])/2
		# print("xl =", x_f[0][0], "xg =", x_f[1][0],"xh =", x_f[2][0])
		# print("f(xl) =", x_f[0][1], "f(xg) =", x_f[1][1],"f(xh) =", x_f[2][1])
		# print("Центр тяжести: хс =", xc)	

		teta = 1
		xNew = project(x_f[2][0], xc, teta)
		l.append(xNew)
		fNew = f(xNew)
		kvcf +=1
		# print("пробное xNew =", xNew)
		# print("пробное f(xNew) =", fNew)
		#compare
		if fNew<x_f[0][1]:
			teta = gamma
		elif fNew<x_f[1][1]:
			teta = alpha
		elif fNew<x_f[2][1]:
			teta = beta[1]
		else:
			teta = beta[0]
			x1 = project(x_f[0][0], x_f[1][0], teta)
			x2 = project(x_f[0][0], x_f[2][0], teta)
			x_f = [(x1,f(x1)), (x2, f(x2)), x_f[0]]
			kvcf +=2
			l.append(x1)
			l.append(x2)
			# print("Редукция")
			i+=1
			continue
		
		xNew = project(x_f[2][0], xc, teta)
		fNew = f(xNew)
		kvcf += 1
		l.append(xNew)
		# Проверка на зацикливание
		# print("prev", previousF[-2], "new", fNew)
		if fNew == previousF[-2] or fNew == previousF[-3] or fNew == previousF[-4]:
			# print("in")
			teta=beta[0]
			xNew = project(x_f[2][0], xc, teta)
			l.append(xNew)
			fNew = f(xNew)
			kvcf += 1
		# print("teta =", teta)
		# print("xNew =", xNew)
		# print("f(xNew) =", fNew, i)
		x_f[2] = (xNew, fNew)
		i+=1
		if i >5000:
			return 1000
		previousF.append(fNew)

	# print([l[i][0] for i in range(len(l))])
	# plt.scatter([l[i][0] for i in range(len(l))], [l[i][1] for i in range(len(l))])
	# plt.scatter(l)
	plt.xlabel('x1')
	plt.ylabel('x2')
	# plt.show()
	
	xc = (x_f[0][0]+x_f[1][0])/2
	# Критерий окончания
	# print("Критерий окончания [1/(n+1)] * sum(f(x)-f(xc))**2 = ",sum([(x_f[i][1]-f(xc))**2 for i in range(n)])/(n+1), "<e")	
	# print("колво редукций:", k)
	x_f.sort(key=lambda p:p[1])	
	# print("xl =", x_f[0][0], "xg =", x_f[1][0],"xh =", x_f[2][0])
	# print("f(xl) =", x_f[0][1], "f(xg) =", x_f[1][1],"f(xh) =", x_f[2][1])
	# print("*хmin =", x_f[0][0])
	# print("*f(хmin) =", x_f[0][1])
	# print("Начальный размер симплекса", simpSizeCoef)
	return i, kvcf

# Функція модифікації методу з випадковими коефіцієнтами
def nelderMidRand(betaRange, gammaRange, simpSizeCoef=0.01):
# Начальные точки
	xStart = array([5, 5])
	x_f = ( xStart-array([0, simpSizeCoef]), None),( xStart+array([simpSizeCoef, 0]), None),( xStart+array([0, simpSizeCoef]) , None)
	x_f = [(p[0], f(p[0])) for p in x_f]
	xc = (x_f[0][0]+x_f[1][0])/2
	startSize = (sum([(x_f[i][1]-f(xc))**2 for i in range(n)])/(n+1))**0.5
	# print("Начальные точки: "+" ".join(["x"+str(i+1)+" = "+str(x_f[i][0]) for i in range(3)]))

	# print(" ".join(["f"+str(i+1)+" = "+str(x_f[i][1]) for i in range(n)]))

	project = lambda xh, xc, teta: xh + (1+teta)*(xc-xh)

	previousF = ['*','*','*', '*', '*']
	k=0
	i=0
	# превычисление (для первой итерации)
	x_f.sort(key=lambda p:p[1])
	while (sum([(x_f[i][1]-f(xc))**2 for i in range(n)])/(n+1))**(0.5) > e:
		# print(str(i+1)+")")
		x_f.sort(key=lambda p:p[1])
		xc = (x_f[0][0]+x_f[1][0])/2
		# print("xl =", x_f[0][0], "xg =", x_f[1][0],"xh =", x_f[2][0])
		# print("f(xl) =", x_f[0][1], "f(xg) =", x_f[1][1],"f(xh) =", x_f[2][1])
		# print("Центр тяжести: хс =", xc)	

		teta = 1
		xNew = project(x_f[2][0], xc, teta)
		fNew = f(xNew)
		# print("пробное xNew =", xNew)
		# print("пробное f(xNew) =", fNew)
		# compare
		if fNew<x_f[0][1]:
			teta = random.choice(gammaRange)
		elif fNew<x_f[1][1]:
			teta = alpha
		elif fNew<x_f[2][1]:
			teta = random.choice(betaRange)
		else:
			teta = -random.choice(betaRange)
			x1 = project(x_f[0][0], x_f[1][0], teta)
			x2 = project(x_f[0][0], x_f[2][0], teta)
			x_f = [(x1,f(x1)), (x2, f(x2)), x_f[0]]
			# print("Редукция")
			i+=1
			k+=1
			continue
		
		xNew = project(x_f[2][0], xc, teta)
		fNew = f(xNew)
		# Проверка на зацикливание
		# print("prev", previousF[-2], "new", fNew)
		if fNew == previousF[-2] or fNew == previousF[-3] or fNew == previousF[-4]:
			# print("in")
			teta=beta[0]
			xNew = project(x_f[2][0], xc, teta)
			fNew = f(xNew)
		# print("teta =", teta)
		# print("xNew =", xNew)
		# print("f(xNew) =", fNew)
		x_f[2] = (xNew, fNew)
		i+=1
		if i >5000:
			return 1000
		previousF.append(fNew)

		
		
	xc = (x_f[0][0]+x_f[1][0])/2
	# Критерий окончания
	# print("Критерий окончания [1/(n+1)] * sum(f(x)-f(xc))**2 = ",sum([(x_f[i][1]-f(xc))**2 for i in range(n)])/(n+1), "<e")	
	# print("колво редукций:", k)
	x_f.sort(key=lambda p:p[1])	
	# print("xl =", x_f[0][0], "xg =", x_f[1][0],"xh =", x_f[2][0])
	# print("f(xl) =", x_f[0][1], "f(xg) =", x_f[1][1],"f(xh) =", x_f[2][1])
	# print("*хmin =", x_f[0][0])
	# print("*f(хmin) =", x_f[0][1])
	# print("Начальный размер симплекса", simpSizeCoef)
	return i	

	
# random mod
gammaRange = np.arange(2, 5, 0.5)
betaRange = np.linspace(0.7, 0.99, len(gammaRange))
# su = 0
# for i in range(100):
	# res = nelderMidRand(betaRange, gammaRange)
	# print(res)
	# su += res
	
# print('Classic result: ', nelderMid(beta, gamma))
# print('Modified result:', su/100) 
	
	

# Дослідження коефіцієнтів
# gamma = np.arange(1.1, 10, 0.2)
# beta = np.linspace(0.05, 0.999, len(gamma))
# alpha = 1
# minres, ming, minb = 1999, 10000, 9
# for b in beta:
	# for g in gamma:
		# res = nelderMid((-b, b),g)
		# if res< minres:
			# minres = res
			# ming, minb = g, b
		# if res==1000:
			# plt.scatter(b,g, 100, marker='x')
		# plt.scatter(b,g,res)
# print(minres, 'if gamma =',ming, ', beta =,', minb)
# plt.show()



# start simplex size
# coefs = np.arange(0.0000001,1, 0.0001)
rLower, rLowerS, rUpper, rUpperS  = 0, 0, 0, 0
min = 1000
# xStart=[5, 5]
sizes = [0.0000001*10**i for i in range(10)]
# sizes = np.arange(0.000001, 1.5, 0.0002)
for c in sizes:
	res = nelderMid1(startSize=c)
	# if res[1]<150:
		# rLower += 1
		# rLowerS += res[1]
	# elif res[1]>150:
		# rUpper +=1
		# rUpperS +=res[1]
	# if res[1]<min:
		# min = res[1]
	print(c, res)
	plt.scatter(c, res[1])

# print('Lower count',rLower, 'Upper count', rUpper)
# lmean = rLowerS/rLower
# umean = rUpperS/rUpper
# print('Lower mean', lmean, 'Upper mean', umean)
# plt.plot([0,3], [lmean, lmean], c='red')
# plt.plot([0,3], [umean, umean], c='red')
# plt.ylabel('КВЦФ')
# plt.xlabel('Початковий розмір симпплексу')
# plt.show()

# Дослідження модифікації з випадковими коефіцієнтами
# for c in coefs:
	# su1 = 0
	# size calc
	# xStart = array([5, 5])
	# x_f = ( xStart-array([0, c]), None),( xStart+array([c, 0]), None),( xStart+array([0, c]) , None)
	# x_f = [(p[0], f(p[0])) for p in x_f]
	# xc = (x_f[0][0]+x_f[1][0])/2
	# startSize = (sum([(x_f[i][1]-f(xc))**2 for i in range(n)])/(n+1))**0.5
	##### nelderMidRand calc
	# for i in range(5):
		# iterations = nelderMidRand(betaRange, gammaRange, c)
		# su1 += iterations
	# plt.scatter(startSize,su1/5, marker='x')
	# plt.scatter(startSize, nelderMid(c), marker='o')
# plt.legend()
# plt.show()



