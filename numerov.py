import numpy as np
import matplotlib.pyplot as plt
def V(x):
	return 0.5*x*x


def numerov(E,plot=0,normalize=0,returnfi=0):
	#lim=5*np.sqrt(2*E)
	xmin = -8.0
	xmax = 8.0
	ndivide = 1000
	
	s=(xmax-xmin)/(ndivide-1)
	
	G = np.zeros(ndivide,dtype=np.dtype('f16'))
	f = np.zeros(ndivide,dtype=np.dtype('f16'))
	x = np.linspace(xmin, xmax, num=ndivide,dtype=np.dtype('f16'))
	
	nn = 0 #number of nodes
	
	#assign initial values of phi
	f[0] = 0
	f[1] = 0.0001
	
	G[0] = 2*V(x[0]) - 2*E
	G[1] = 2*V(x[1]) - 2*E
	
	for i in range(2, ndivide):
		G[i] = 2*V(x[i]) - 2*E
		f[i] = (-f[i-2] + 2*f[i-1] + 5.0*G[i-1]*f[i-1]*s*s/6.0 + G[i-2]*f[i-2]*s*s/12.0)/(1-G[i]*s*s/12.0)
		if (f[i]*f[i-1]) < 0.0:
			nn+=1
	
	#print('E = %.10f Phimax = %.5f nn = %d'%(E,f[i],nn))
	if normalize:
		I=np.trapz(x=x,y=f**2)
		f=f/np.sqrt(I)
	if plot:
		plt.plot(x,f,label='%.10f'%E)
		#plt.savefig('%.10f.png'%E)
		#plt.clf()
	if returnfi:
		return nn,f[i]
	else:
		return nn

#uniform sampling, found limits of nn changes
Es=np.linspace(0,7,num=200,dtype=np.dtype('f16'))
NNs=[numerov(i) for i in Es]
#print(NNs)

Pairs=[]
for i in range(1,len(NNs)):
	if NNs[i]!=NNs[i-1]:
		Pairs.append((Es[i-1],Es[i],NNs[i-1],NNs[i]))

epsilon=0.0001
print('With accuracy of %.10f'%epsilon)
for a,b,nn0,nn1 in Pairs:
	#print(a,b,nn0,nn1)
	x0=a
	x1=b
	_,fxnew=numerov(x0,normalize=1,returnfi=1)
	while(abs(fxnew) > epsilon):
		xnew=0.5*(x0+x1)
		nn_new,fxnew=numerov(xnew,normalize=1,returnfi=1)
		print(xnew,nn_new,fxnew)
		if (nn_new == nn0):
			x0 = xnew
		elif (nn_new == nn1):
			x1 = xnew
	print('Found eigenenergy at %.10f, nu = %d'%(x0,nn0))
	numerov(x0,plot=1,normalize=1)

plt.legend()
plt.xlabel(r'$x_r$')
plt.title('First several wavefunctions')
plt.savefig('sum.png')
