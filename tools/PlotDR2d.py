from pylab import *
from sys import argv

dat=genfromtxt(argv[1])
n=len(dat)
sn=int(sqrt(n))
w=dat[:,0].reshape(sn,sn)
g=dat[:,1].reshape(sn,sn)
DR2d=dat[:,2].reshape(sn,sn)

contourf(w,g,log10(DR2d),20,cmap='gist_heat')
xlabel('Frequency')
ylabel('Damping rate')
colorbar()
tight_layout()
show()
