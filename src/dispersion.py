# Author: Jonathan Cookmeyer, Daniel Told
# Summer 2015

from __future__ import division,absolute_import,print_function,unicode_literals #for Python 2.7

#from pylab import *
from scipy.special import wofz,jn,ive
from scipy.optimize import fsolve,root,newton

from numpy import amin,matrix,tile,linspace,repeat,empty,log10,sqrt,exp,pi
from numpy.linalg import det,eig
#from matplotlib.pyplot import contourf,colorbar,show
#import matplotlib.pyplot as plt
import numpy as np
from sys import exit

#Define plasma dispersion function (scaling factor 1j*sqrt(pi) is multiplied separately
#for performance reasons)
def Z(zeta):
    return wofz(zeta)

#Initialize some variables only once
def init(params,outfile_obj):
    global first, first_info, sqrtpi1j, outfile, warned_Bessel
    #make outfile global for this module
    outfile=outfile_obj
    first=True
    first_info=[]

    sqrtpi1j=1j*sqrt(pi)
    warned_Bessel=0.

def init_params(params,silent=False):
    global kperp,kpar,beta,tau,Tpar_Tperp,gam,eta,nb,theta,k,first
    first=True
    wavevector_mode=params[0]
    kperp=params[1]
    kpar=params[2]
    beta=params[3]
    tau=params[4]
    Tpar_Tperp=params[5]
    gam=params[6] #1 or 5/3
    eta=params[7]
    # discretization parameters
    nb=params[8] #how far to take sum
    theta=params[9]
    k=params[10]
    if wavevector_mode==2:
        kpar=k*np.cos(theta*pi/180)
        kperp=k*np.sin(theta*pi/180)
        if not silent: print('Running with k_par=%7.4g and k_perp=%7.4g.' %(kpar,kperp),file=outfile)


#Generate Bessel function array -- this is only done once, unless kperp or Tpar_Tperp change
def make_J_array(kperp,v,nb,Tpar_Tperp):
    global jn_arr, jnp1_arr, jnp2_arr
    jn_arr = np.empty((2*nb+3),dtype=np.float64)
    kv_t=kperp*v/sqrt(Tpar_Tperp)
    
    js=np.arange(0,2*nb+3)
    j_shift=js-nb
    for j in js:
        jn_arr[j] = jn(j_shift[j],kv_t)

#Generate Bessel function array -- this is only done once, unless kperp or Tpar_Tperp change
def make_Ive_array(kperp,nb,Tpar_Tperp):
    global ive_arr
    ive_arr = np.empty((2*nb+2),dtype=np.float64)
    b=0.5*kperp**2/Tpar_Tperp
    
    js=np.arange(0,2*nb+2)
    j_shift=js-nb-1
    for j in js:
        ive_arr[j] = ive(j_shift[j],b)

#Generate plasma dispersion function array -- this has to be done any time omega or kpar
#changes, i.e. also when the root finder moves
def make_Z_array(omega,kpar,nb):
    global zn_arr, znp1_arr, znm1_arr

    arg_Z=omega-np.arange(-nb-1,nb+2)/kpar
    zn_arr=Z(arg_Z)*sqrtpi1j

#This routine manages the generation of the Bessel and Z function arrays for analytical vperp integrals
def calc_bessel_Z_analytical(j, omega, kpar, kperp, Tpar_Tperp, nb):
    global first, first_info
    if first:
        make_Z_array(omega,kpar,nb)
        make_Ive_array(kperp,nb,Tpar_Tperp)
        first=False
        first_info=[kpar,kperp,Tpar_Tperp,omega]
    if first_info[0]!=kpar or first_info[3]!=omega:
        make_Z_array(omega,kpar,nb)
        first=False
        first_info=[kpar,kperp,Tpar_Tperp,omega]
    if first_info[1]!=kperp or first_info[2]!=Tpar_Tperp:
        make_Ive_array(kperp,nb,Tpar_Tperp)
        first=False
        first_info=[kpar,kperp,Tpar_Tperp,omega]
    j_shift=j+nb+1
    return zn_arr[j_shift+1],zn_arr[j_shift-1],zn_arr[j_shift],ive_arr[j_shift],ive_arr[j_shift-1]

#Compute the term with number "j" for all Bessel sums
def calc_sums_analytical(j, omega, kpar, kperp, znp1,znm1,zn,iveval,ivevalm1):
    b=0.5*kperp**2/Tpar_Tperp
#Define some abbreviations to reduce number of necessary operations
    inv_b=1./b
    inv_kpar=1./kpar
    gn=iveval 
    gnm1=ivevalm1
    gnpr=gnm1-(j*inv_b+1)*gn
    bgnpr1=b*gnm1
    bgnpr2=-gn*(j+b)
    bgnpr3=-gn*b
    gn_aux=gn+bgnpr1+bgnpr2

    omj_zn=(omega-float(j)*inv_kpar)*zn
    omjp1_znp1=(omega-(1.+j)*inv_kpar)*znp1
    omjm1_znm1=(omega-(-1.+j)*inv_kpar)*znm1

    summand=np.zeros(13,np.complex128)
    summand[0]=znm1*(bgnpr1+bgnpr3)                       #the m-1 half of I_1
    summand[1]=2*j*znp1*gn-znp1*(bgnpr1+bgnpr3)           #the m+1 half of I_1
    summand[2]=(1+omjm1_znm1)*(bgnpr1+bgnpr3)
    summand[3]=2*j*(1+omjp1_znp1)*gn-(1+omjp1_znp1)*(bgnpr1+bgnpr3)
    summand[4]=(1+omj_zn)*gn
    summand[5]=znm1*gn_aux
    summand[6]=znp1*gn_aux
    summand[7]=(1+omjm1_znm1)*gn_aux
    summand[8]=(1+omjp1_znp1)*gn_aux
    summand[9]=(j-1)*znm1*gnpr
    summand[10]=j*(j-1)*znm1*gn
    summand[11]=(j-1)*(1+omjm1_znm1)*gnpr
    summand[12]=j*(j-1)*(1+omjm1_znm1)*gn

    return summand

def get_sums_analytical(kperp, kpar, omega, Tpar_Tperp, nb):
    global warned_Bessel
    #minimum accuracy 
    acc=1e-12

    #zeroth term
    j=0
    znp1,znm1,zn,iveval,ivevalm1=calc_bessel_Z_analytical(j, omega, kpar, kperp, Tpar_Tperp, nb)
    totalsum=calc_sums_analytical(j, omega, kpar, kperp, znp1, znm1, zn,iveval,ivevalm1)
    for j in range(1,nb+1):
        #positive j
        znp1,znm1,zn,iveval,ivevalm1=calc_bessel_Z_analytical(j, omega, kpar, kperp, Tpar_Tperp, nb)
        psummand=calc_sums_analytical(j, omega, kpar, kperp, znp1, znm1, zn,iveval,ivevalm1)
        totalsum+=psummand

        #negative j
        znp1,znm1,zn,iveval,ivevalm1=calc_bessel_Z_analytical(-j, omega, kpar, kperp, Tpar_Tperp, nb)
        nsummand=calc_sums_analytical(-j, omega, kpar, kperp, znp1, znm1, zn,iveval,ivevalm1)
        totalsum+=nsummand

        if any(np.isnan(totalsum)): 
            print('WARNING: Encountered NaN in Bessel sums! Entries: ',j,omega,np.isnan(totalsum),zn,znp1,znm1,file=outfile)
            exit()
        #if necessary accuracy has been reached, truncate at present i
        if all(abs(psummand+nsummand)/abs(totalsum)<acc) and j>2:
#            print ('Used %d terms in Bessel sums.' %j,file=outfile)
            break
        else:
            if j==nb-1 and warned_Bessel!=kperp:
                print ('Warning: Could not reach desired convergence accuracy in Bessel sums! Max. deviation: %g' %(max(abs(psummand+nsummand)/abs(totalsum))),file=outfile)
                warned_Bessel=kperp
    return list(totalsum)

#new version computation using analytically evaluated vperp integrals
def dispersion_relation_analytical(omega,beta,tau,Tpar_Tperp,kperp,kpar,gam,eta,nb,theta,k):
    k2=kperp**2+kpar**2
    b=0.5*kperp**2/Tpar_Tperp
    inv_kpar=1./kpar
    inv_kperp=1./kperp
    inv_b=1./b
    summand=get_sums_analytical(kperp,kpar,omega,Tpar_Tperp,nb)
    M = 1j*eta*k2*inv_kpar + omega - 0.5*inv_kpar*(1j*eta*kperp**2*inv_kpar+omega)*(summand[6]-summand[5]) + 0.5*(Tpar_Tperp-1.)/Tpar_Tperp*inv_kpar*(summand[8]-summand[7]) + 1j*eta*inv_kpar*(summand[2]-summand[3])
    N = 1j/beta*k2*inv_kperp + 1j*omega*inv_kperp*(summand[9]+inv_b*summand[10]-0.5*(summand[6]+3*summand[5])) - 1j*inv_kperp*(Tpar_Tperp-1.)/Tpar_Tperp*(summand[11]+inv_b*summand[12]-0.5*(summand[8]+3*summand[7]))
    O = 0.5j*gam/tau*(-inv_kperp*(summand[2]-summand[3]) + 0.5*kperp*inv_kpar*(summand[6]-summand[5]))
    
    P = 1j*kpar/beta - 1j*inv_kpar*(1j*eta*kperp**2*inv_kpar+omega)*(summand[9]+inv_b*summand[10]+0.5*(summand[6]-summand[5])) + 1j*inv_kpar*(Tpar_Tperp-1.)/Tpar_Tperp*(summand[11]+inv_b*summand[12]+0.5*(summand[8]-summand[7])) + eta*inv_kpar*(summand[2]+summand[3])
    Q = -(1j*eta*k2+omega*kpar)*inv_kperp + 0.5*omega*inv_kperp*(summand[6]-summand[5]) - 0.5*(Tpar_Tperp-1.)/Tpar_Tperp*inv_kperp*(summand[8]-summand[7])
    R = -0.5*gam/tau*(inv_kperp*(summand[2]+summand[3]) + kperp*inv_kpar*(summand[9] + inv_b*summand[10] + 0.5*(summand[6]-summand[5])))

    S = -(1j*eta*kperp**2*inv_kpar+omega)*1j*inv_kperp*inv_kpar*Tpar_Tperp*(summand[0]+summand[1]) + 1j*inv_kpar*inv_kperp*(Tpar_Tperp-1)*(summand[2]+summand[3]) + 2*eta*kperp*inv_kpar*summand[4]
    T = -Tpar_Tperp*omega*inv_kperp**2*(summand[0]-summand[1]) + inv_kperp**2*(Tpar_Tperp-1)*(summand[2]-summand[3])
    U = -1 - 0.5*gam/tau*(2*summand[4]+Tpar_Tperp*inv_kpar*(summand[0]+summand[1]))

    global mat
    mat=[[M,N,O],[P,Q,R],[S,T,U]]
    if det(mat)>1: print(summand[6],-summand[5],file=outfile)
    return det(mat)

             
#wrapper for analytical or numerical dispersion relation
def DR(*arg):
    return dispersion_relation_analytical(*arg)

#dispersion relation evaluation used by the root finding routine
def DR_RF(x):
        #For very negative imaginary parts, sometimes NaNs in the plasma dispersion function
        #occur. We just return 0,0 for gamma<-20, stopping the root solver. If all goes well,
        #the outer routines should discover the discontinuity of that solution and try again 
        #at a different point, hopefully circumventing this problem.
        if x[1]<-20:
            return [1e9,1e9]
        omega=x[0]+1j*x[1]
        out=DR(omega,beta,tau,Tpar_Tperp,kperp,kpar,gam,eta,nb,theta,k)
        return [out.real,out.imag]
def DR_1d(x):
	omega=x
	out=DR(omega,beta,tau,Tpar_Tperp,kperp,kpar,gam,eta,nb,theta,k)
	return out.real

#Compute dispersion relation for a single point and return the determinant
def DR_point(params,x):
        init_params(params,True)
        omega=x[0]+1j*x[1]
        return abs(DR(omega,beta,tau,Tpar_Tperp,kperp,kpar,gam,eta,nb,theta,k))

#Runs the root finder from a given start value and returns the result
def DR_Solve(params,start):
        init_params(params)

        solt=root(DR_RF,x0=start,method='lm',tol=1e-12)
        if solt.get('success')==False:
            x, infodict, ier, mesg = fsolve(DR_RF, x0=start,maxfev=0,
				            full_output=1,xtol=1e-6, epsfcn=1e-6)
            nfev=infodict['nfev']
        else:
            x = solt.get('x')
            ier=1
            mesg=solt.get('message')
            nfev=solt.get('nfev')
        if ier:
            print('Converged after', nfev, 'iterations',file=outfile)
            print('\n RESULT:',(x[0],x[1]),'\n',file=outfile)
            residual=DR(x[0]+1j*x[1],beta,tau,Tpar_Tperp,kperp,kpar,gam,eta,nb,theta,k)
            print('DR residual:', residual,file=outfile)
            if abs(residual)>1e-8: 
                print('WARNING: Large residual, root may not be properly converged!',file=outfile)
                return None, None, None, None, 0,residual
                #convergence_issue+=1
            omega_z = x[0]+1j*x[1]
            compeigvec=None #for complex components
            u, v = eig(mat)
            eps=2*amin(abs(u))
            move=abs(x[0]-start[0])+abs(x[1]-start[1]) #will be 0 if it doesn't move

	    #
	    # Get the eigenvectors.
	    #
            print('Smallest eigenvalue',eps,file=outfile)
            By=None
            dn=None
            for i in range(3):
                if abs(u[i]) < eps and eps < 1e-4:
            	    #print ("Smallest eigenvalue and corresponding eigenvector:",file=outfile)
            	    #print (abs(u[i]),file=outfile)
                    compeigvec=np.array([0+1j*0,0+1j*0,0+1j*0])
                    for j in range(3):
                        compeigvec[j]=v[j,i]
                    By=compeigvec[0]/compeigvec[1]
                    dn=compeigvec[2]/compeigvec[1]
            return x[0],x[1],By,dn,1,residual
        else:
            return None, None, None, None, 0,-1

def PlotDR2d(params,np=(100,100),cntr=[0,-1],raxis=3.,iaxis=3.,outfilename='output.log'):
    global kperp,beta,nb,theta,first
    outfile=open(outfilename,'w')
    init(params,outfile)
    init_params(params)

    xmin = cntr[0] - raxis/2; xmax = cntr[0] + raxis/2
    ymin = cntr[1] - iaxis/2; ymax = cntr[1] + iaxis/2
    arr1 = tile(linspace(xmin, xmax, np[0]), np[1])
    arr2 = repeat(linspace(ymin, ymax, np[1]), np[0])
    omega_out = (arr1 + 1j * arr2).reshape(np[1], np[0])
    omega = (arr1 + 1j * arr2)

    
    DR2d=empty((np[0]*np[1]),dtype=float)
    for i in range((np[0]*np[1])):
        if (i % np[0])==0:
            print('Computing line ', int(i/np[0])+1, 'of', np[1],file=outfile)
            outfile.flush()
        DR2d[i] = abs(DR(omega[i],beta,tau,Tpar_Tperp,kperp,kpar,gam,eta,nb,theta,k))
    DRout = DR2d.reshape(np[1], np[0])

    return omega_out,DRout


