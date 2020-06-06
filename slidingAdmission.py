import matplotlib.pyplot as plt
from numpy import array
import random
from math import sqrt
import numpy as np
from nelderMid import nelderMid, nelderMid1

def areEqual(arr1, arr2): 
  
    # If lengths of array are not  
    # equal means array are not equal 
	n, m = len(arr1), len(arr2)
	if (n != m): 
		return False; 
	for i in range(0, n - 1): 
		if (arr1[i] != arr2[i]): 
			return False; 
	return True; 

# Функция Розенброка
# a, b = 1, 100
# f = lambda x: (a-x[0])**2 + b*(x[1]-x[0]**2)**2

# f = lambda x: 4*x[0]-x[1]**2-12
# g = [ lambda x: x[0]+2, lambda x: -x[0]+0.5, lambda x: x[1]+3]
# h = []
f1 = lambda x: -2
f2 = lambda x: 2
f3 = lambda x: 0.5
# f3 = lambda x: x**2
# g = [lambda x: x[1]-2 ]
# h = [lambda x: 25-x[0]**2-x[1]**2]
# g = [lambda x: 10*x[0]-x[0]**2+10*x[1]-x[1]**2-34, lambda x: x[1], lambda x: x[0]]
f = lambda x: x[0]**2 + x[1]**2
# f1 = lambda x: sqrt(16-(x-4.5)**2)
# f1 = lambda x: 0
# f2 = lambda x: sqrt(25-x**2)
# f3 = lambda x: -sqrt(16-(x-5)**2)+5
# Начальные парараметры
g=[]
h = [lambda x: x[0]**2 + x[1]**2 - 9* x[1] + 4.25]
n = 2
t = 1
m = len(h)
r = n - m
alpha, beta, gamma = 1, (-0.5, 0.5), 2
e = 10e-8

Phi0 = lambda m, t:2*(m+1)*t
Phi = lambda phiPrev, m, r, x, xc: min(phiPrev, (m+1)/(r+1) * sum([sum((xi-xc)**2) for xi in x]))
T = lambda x, h=h, g=g: abs( (sum([hi(x)**2 for hi in h]) + sum([gi(x)**2 if gi(x) < 0 else 0 for gi in g]))**(1/2) )



xStart = (4, 4.5)
def slidingDecent(e=e):
	kvcf = 0
	x_f = (array([3.592, 4.092]), None),(array([4.558, 4.351]), None),(array([3.85, 5.06]) , None)
	# x_f = (array([-2, 3]), None),(array([-2.25, 3.25]), None),(array([-1.75, 3.75]) , None)
	# x_f = (array([1, 1]), None),(array([1.057, 1.015]), None),(array([1.015, 1.057]) , None) # = + >=
	# x_f = (array([-1.9, 2.5]), None),(array([-1.8, 2.6]), None),(array([-1.85, 2.44]) , None)
	x_f = [(p[0], f(p[0])) for p in x_f]
	kvcf += 3
	xTrack = [x_f[i][0] for i in range(len(x_f))]
	i = 0
	phi = (Phi0(m, t))
	phiPrev=-1
	# 
	while phi>e and not areEqual(x_f[1][0], x_f[0][0]):
	# while kvcf<600:
		x_f.sort(key=lambda p:p[1])
		xc = (x_f[0][0] + x_f[1][0])/2
		tt = T(x_f[0][0], h, g)
		if phi-tt > 0:
			pass
		else:
			# min T(X) start from xNew
			xNew = x_f[0][0]
			xMinforT, incKvcf = nelderMid(phi, f=T, xStart=xNew)
			x_f[2] = (xMinforT, f(xMinforT))
			x_f.sort(key=lambda p:p[1])
			xc = (x_f[0][0] + x_f[1][0])/2
			xTrack = xTrack + [x_f[2][0]]
			kvcf +=incKvcf + 1  
		teta = 1
		project = lambda xh, xc, teta: xh + (1+teta)*(xc-xh)
		# x test
		xNew = project(x_f[2][0], xc, teta)
		fNew = f(xNew)
		kvcf += 1
		tt = T(xNew, h, g)
		if phi-tt > 0:
			pass
		else:
			# min T(X) start from xNew
			xNew, incKvcf = nelderMid(phi, f=T, xStart=xNew)
			fNew = f(xNew)
			kvcf += incKvcf
			xTrack.append(xNew)
		x_f.sort(key=lambda p:p[1])
		#compare
		if fNew<x_f[0][1]:
			teta = gamma
		elif fNew<x_f[1][1]:
			teta = alpha
		elif fNew<x_f[2][1]:
			teta = beta[1]
		else:
			teta = beta[0]
		if teta==beta[0]:
			x1 = project(x_f[0][0], x_f[1][0], teta)
			x2 = project(x_f[0][0], x_f[2][0], teta)
			x_f = [(x1,f(x1)), (x2, f(x2)), x_f[0]]
			kvcf += 2
			x_f.sort(key=lambda p:p[1])
			
			xNew = x_f[0][0]
			xTrack.append(x1)
			xTrack.append(x2)
			thing = [x_f[i][0] for i in range(len(x_f))] 
			xc = sum(thing)/len(x_f)
			phi = Phi(phi, m, r, thing, xc)
			i+=1
			continue
		else:
			xc = (x_f[0][0] + x_f[1][0])/2
			xNew = project(x_f[2][0], xc, teta)
		tt = T(xNew, h, g)
		if phi-tt > 0:
			x_f[2] = (xNew, f(xNew))
			xTrack.append(xNew)
			kvcf += 1
			x_f.sort(key=lambda p:p[1])
		else:
			# min T(X) start from xNew
			xMinforT, incKvcf = nelderMid(phi, f=T, xStart=xNew)
			x_f[2] = (xMinforT, f(xMinforT))
			kvcf += incKvcf + 1
			xTrack.append(xMinforT)
			x_f.sort(key=lambda p:p[1])
		if areEqual(xTrack[-1], xTrack[-2]) or areEqual(xTrack[-1], xTrack[-3]):
			teta=beta[1]
			x1 = project(x_f[0][0], x_f[1][0], teta)
			x2 = project(x_f[0][0], x_f[2][0], teta)
			x_f = [(x1,f(x1)), (x2, f(x2)), x_f[0]]
			kvcf += 2
			x_f.sort(key=lambda p:p[1])
			
			xNew = x_f[0][0]
			xTrack.append(x1)
			xTrack.append(x2)
			thing = [x_f[i][0] for i in range(len(x_f))] 
			xc = sum(thing)/len(x_f)
			phi = Phi(phi, m, r, thing, xc)
			i+=1
			continue
		i+=1
		thing = [x_f[i][0] for i in range(len(x_f))]
		if i > 1000:
			break
		xc = sum(thing)/len(x_f)
		phiPrev = phi
		phi = Phi(phi, m, r, thing, xc)
		print(phi)
	
	x = [xTrack[i][0] for i in range(len(xTrack))]
	y = [xTrack[i][1] for i in range(len(xTrack))]
	plt.scatter(x,y)
	ar = np.arange(-2,3,0.01)
	plt.plot(ar, [f1(ari) for ari in ar], c='yellow')
	plt.plot([f1(ari) for ari in ar], ar, c='yellow')
	ar = np.arange(-2,3,0.01)
	plt.plot([f3(ari) for ari in ar], ar, c='yellow')
	ar = np.arange(-2,2,0.01)
	plt.plot(ar, [f2(ari)+1 for ari in ar], c='yellow')
	print('i =', i)
	plt.xlabel('x1')
	plt.ylabel('x2')
	plt.show()
	
	return phi, kvcf
	
	
	
	
for i in range(1,10):
	print(slidingDecent(e=0.1**i))











