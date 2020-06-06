def calc(x):
	restartVal = 90
	if x[-1][1]<0.00001:
		return []
	return x[:-1] + calc([(restartVal+x[0][0], x[-1][1]*0.44), (restartVal+x[-1][0], 0.56*x[1][1])])
	
	
a = [(70.185,0.44), (248.812, 0.56)]
print(a[0])
print(a[1])
res = calc(a)
s, Mnew = 0, 0
for t in res:
	print(t)
	s+=t[1]
	Mnew += t[0]*t[1]
print('cумма', s)
print('Класичний: М =', a[0][0]*a[0][1]+a[1][0]*a[1][1] )
print('Модифікований: М =',Mnew)