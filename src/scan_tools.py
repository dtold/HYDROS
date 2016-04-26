from __future__ import division, print_function

#import matplotlib.pyplot as plt
import numpy as np
from numpy import log10
from sys import  stderr, stdout, exit
from dispersion import DR_Solve, init, DR_point, PlotDR2d
from scipy.interpolate import InterpolatedUnivariateSpline as US

def banner():
        line='_'*80+'\n'
        print(line)
        print(' '*30+'THIS IS HYDROS v1.0.0')
        print(line)

# This function allows the roots to be found ranging across XARRAY
def HYDROS(params,start,scan1,scan2,log=1,method='default',scan_options=[1000.,1],
                 outfilename='output.log',scanfilename='scan.log',scan_type='1d',normalization=1):
        banner()
        global outfile
        check_variables(params,scan1,scan2)
        
        outfile=open(outfilename,'w')
        scanfile=init_resultfile(scanfilename)

        if scan_type=='DR-map':
                w,DR2d=PlotDR2d(params,np=(800,800),cntr=[1.5,-1.5],raxis=3.,iaxis=3.,outfilename='output.log')
                write_2d_scan('DR2d.log',w.real,w.imag,DR2d)
                return
        
        scanvar,precis,scan_range,log,scan_list=scan1
        scansize=precis

        scanvar2,precis2,scan_range2,log2,scan_list2=scan2
        if scan_type=='2d':
                scansize*=precis2
        else:
                if scan_list2==[]:
                        scanvar2,precis2,scan_range2,log2,scan_list2=[None,None,None,None,[]]
        
        #initialize dispersion relation solver
        params,start,scan_range,scan_range2=renormalize(normalization,scanvar,scanvar2,params,start,scan_range,scan_range2,scan_type)
        init(params,outfile)

        #set up scan arrays
        XARRAY, X2ARRAY=setup_xarray(scan_list,scan_list2,scan_range,scan_range2,log,log2,precis,precis2,scan_type)

        OMEGA=np.zeros(scansize)
        GAMMA=np.zeros(scansize)
        BY=np.zeros(scansize,np.complex)
        DN=np.zeros(scansize,np.complex)
        i=0
        while i<len(XARRAY): 
                outfile.flush()
                update_output('\r Computing index %d/%d.' %(i+1,len(XARRAY)))
                print("\n=====================================================================",file=outfile) 
                if X2ARRAY!=None:
                        print("\nComputing index no. %d, %s = %7.3g, %s = %7.3g" %(i+1, scanvar, XARRAY[i], scanvar2, X2ARRAY[i]),file=outfile) 
                else:
                        print("\nComputing index no. %d, %s = %7.3g" %(i+1, scanvar, XARRAY[i]),file=outfile) 

                if scan_type=='2d' and i%precis==0 and i>0: 
                        start=[OMEGA[i-precis],GAMMA[i-precis]]
                        shift=precis-1
                else:
                        shift=0

                scan=call_DR_get_root(i,params,start,scanvar,scanvar2,method,scan_options,XARRAY,OMEGA,GAMMA,BY,DN,X2ARRAY,shift)

                omega,gamma,By,dn,err=scan

                if err==1: # no error
                        OMEGA[i]=omega
                        GAMMA[i]=gamma
                        BY[i]=By
                        DN[i]=dn
                else:
                        update_output('\r Lost root on index %d/%d, aborting scan.' %(i+1,len(XARRAY)),True)
                        print("Error: Lost root, ending scan at scanvar=",XARRAY[i],file=outfile)
                        if scan_type=='2d':
                                return XARRAY[:i], OMEGA[:i], GAMMA[:i], BY[:i], DN[:i], X2ARRAY[:i]
                        else:
                                return XARRAY[:i], OMEGA[:i], GAMMA[:i], BY[:i], DN[:i], None
                if X2ARRAY!=None:
                        b=X2ARRAY[i]
                else:
                        b=0.
                a,b,c,d=renormalize_output(normalization,scanvar,scanvar2,params,XARRAY[i],b,OMEGA[i],GAMMA[i])
                write_result_to_scanfile(scanfile,a,b,c,d,BY[i],DN[i])
                i+=1
        update_output('\r Scan complete. Logfile: \'%s\', Scan output: \'%s\'.' %(outfilename,scanfilename),True)
        scanfile.close()

def check_variables(params,scan1,scan2):
        scanvar,precis,scan_range,log,scan_list=scan1 
        scanvar2,precis2,scan_range2,log2,scan_list2=scan2
        wavevector_mode=params[0]
        kperp=params[1]
        kpar=params[2]
        theta=params[9]
        k=params[10]
        if wavevector_mode==1:
                if scanvar in ['k','theta']:
                        exit('Error: scanning over k or theta in wavevector_mode = 1 makes no sense!')
        if wavevector_mode==2:
                if scanvar in ['kpar','kperp']:
                        exit('Error: scanning over kpar or kperp in wavevector_mode = 2 makes no sense!')
     
    
def get_var_from_scan_range(var,value,scanvar,scanvar2,scan_range,scan_range2,scan_type):
        if scanvar==var:
                var=scan_range[0]
        if scan_type=='2d' and scanvar2==var:
                var=scan_range[0]
        return value
                

#reset some parameters to match internal normalization of HYDROS (k_norm=k*rho_i||, omega_norm=omega/kpar/v_Ti||)
def renormalize(normalization,scanvar,scanvar2,params,start,scan_range,scan_range2,scan_type):
        wavevector_mode=params[0]
        kperp=params[1]
        kpar=params[2]
        theta=params[9]
        k=params[10]
        if normalization==1:
                if scanvar=='beta' or scanvar2=='beta':
                        exit('Error: d_i normalization presently not possible for beta scans.')
                beta=params[3]
                #we need kpar to change the normalization of the start frequency
                #if scans over k and/or theta are done, have to overwrite the above values
                if wavevector_mode==1:
                        kpar=get_var_from_scan_range('kpar',kpar,scanvar,scanvar2,scan_range,scan_range2,scan_type)
                        kperp=get_var_from_scan_range('kperp',kperp,scanvar,scanvar2,scan_range,scan_range2,scan_type)
                if wavevector_mode==2:
                        k=get_var_from_scan_range('k',k,scanvar,scanvar2,scan_range,scan_range2,scan_type)
                        theta=get_var_from_scan_range('theta',theta,scanvar,scanvar2,scan_range,scan_range2,scan_type)
                if wavevector_mode==1 and (kperp<0. or kpar<0): exit('Need positive kpar and kperp for wavevector_mode=1!') 
                if wavevector_mode==2:
                        if (k<0. or theta<0): exit('Need positive k and theta for wavevector_mode=2!')
                        elif theta>0. and k>0:
                                kpar=k*np.cos(theta*np.pi/180)
                #k is now given in terms of d_i units, convert it to code internal rho_i|| units
                kpar*=np.sqrt(beta)
                #change start frequency normalization
                start=list(np.array(start,dtype=float)/kpar)

                #now redefine scan ranges to code internal units
                if scanvar in ['kpar','kperp','k']:
                        scan_range=tuple(np.array(scan_range,dtype=float)*np.sqrt(beta))
                if scanvar2 in ['kpar','kperp','k']:
                        scan_range2=tuple(np.array(scan_range2,dtype=float)*np.sqrt(beta))
                #finally, adapt global parameters
                params[1]*=np.sqrt(beta)
                params[2]*=np.sqrt(beta)
                params[10]*=np.sqrt(beta)
        return params,start,scan_range,scan_range2

#switch from internal normalization to output normalization
def renormalize_output(normalization,scanvar,scanvar2,params,a,b,c,d):
        wavevector_mode=params[0]
        beta=params[3]
        theta=params[9]
        kpar=params[2]
        if scanvar=='kpar':
                kpar=a
        if scanvar2=='kpar':
                kpar=b
        if scanvar=='theta':
                theta=a
        if scanvar2=='theta':
                theta=b
        if scanvar=='k':
                k=a
        if scanvar2=='k':
                k=b
        if wavevector_mode==2 and theta>0 and k>0: kpar=k*np.cos(theta*np.pi/180)
        
        if normalization==1:
                c*=kpar
                d*=kpar
                if scanvar in ['kpar','kperp','k']:
                        a/=np.sqrt(beta)
                if scanvar2 in ['kpar','kperp','k']:
                        b/=np.sqrt(beta)
        return a,b,c,d
                        
def init_resultfile(resultfilepath):
        file=open(resultfilepath,'w')
        file.write('#%15s %16s %16s %16s %16s %16s %16s %16s\n'% ('var1','var2','omega','-gamma','Re(dBy/dBz)','Im(dBy/dBz)','Re(dn/dBz)','Im(dn/dBz)'))
        return file

def write_result_to_scanfile(file,a,b,c,d,e,f):
        file.write('%16.8e %16.8e %16.8e %16.8e %16.8e %16.8e %16.8e %16.8e\n'%(a,b,c,d,e.real,e.imag,f.real,f.imag))

def write_2d_scan(file,w,g,DR):
        outfile=open(file,'w')
        outfile.write('#%15s %16s %16s\n' %('omega','gamma','DR'))
        for i in range(len(w.flatten())):
                outfile.write('%16.8e %16.8e %16.8e\n' %(w.flatten()[i],g.flatten()[i],DR.flatten()[i]))
        outfile.close()
                
def setup_xarray(scan_list,scan_list2,scan_range,scan_range2,log,log2,precis,precis2,scan_type):
        if scan_type=='1d':
                #for scanning a predefined list of numbers
                if scan_list!=[]:
                        #interpolate to prescribed length
                        l1_spl=US(range(len(scan_list)),scan_list)
                        XARRAY=l1_spl(np.linspace(0,len(scan_list)-1,precis))
                        X2ARRAY=None

                        #for having a second simultaneously varying variable
                        if scan_list2!=[]:
                                #interpolate to prescribed length
                                l2_spl=US(scan_list,scan_list2)
                                X2ARRAY=l2_spl(XARRAY)

                #log- or linearly spaced 1d scan
                if (scan_list==[] and scan_list2==[]):
                        X2ARRAY=None #np.empty(1,dtype=np.float)
                        if log:
                                XARRAY=np.logspace(np.log10(scan_range[0]),np.log10(scan_range[1]),precis)
                        else:
                                XARRAY=np.linspace(scan_range[0],scan_range[1],precis)
        elif scan_type=='2d':
                #for scanning a predefined list of numbers
                if scan_list!=[]:
                        #interpolate to prescribed length
                        l1_spl=US(range(len(scan_list)),scan_list)
                        XARRAY=np.tile(l1_spl(np.linspace(0,len(scan_list)-1,precis)),precis).flatten()

                        #for having a second simultaneously varying variable
                        if scan_list2!=[]:
                                #interpolate to prescribed length
                                l2_spl=US(scan_list,scan_list2)
                                X2ARRAY=np.repeat(l2_spl(XARRAY),precis)
                        else:
                                exit('Error: need to specify two scan_lists for 2d scan!')

                #log- or linearly spaced 2d scan
                if (scan_list==[] and scan_list2==[]):
                        if log:
                                XARRAY=np.tile(np.logspace(np.log10(scan_range[0]),np.log10(scan_range[1]),precis),precis2).flatten()
                        else:
                                XARRAY=np.tile(np.linspace(scan_range[0],scan_range[1],precis),precis2).flatten()
                        if log2:
                                X2ARRAY=np.repeat(np.logspace(np.log10(scan_range2[0]),np.log10(scan_range2[1]),precis2),precis).flatten()
                        else:
                                X2ARRAY=np.repeat(np.linspace(scan_range2[0],scan_range2[1],precis2),precis).flatten()
        return XARRAY, X2ARRAY

def update_output(string,end=False):
        stdout.write(string.ljust(80))
        stdout.flush()
        if end: stdout.write('\n\n')


def call_DR_get_root(i,params,start,scanvar,scanvar2,method,scan_options,XARRAY,OMEGA,GAMMA,BY,DN,X2ARRAY,shift):
	#no gam scans because only takes value of 1 or 5/3
        vardict={'kperp': 1, 'kpar': 2, 'beta': 3, 'tau': 4, 'Tpar_Tperp': 5, 'eta': 7,
                 'nb': 8, 'theta': 9, 'k': 10}
        params[vardict[scanvar]]=XARRAY[i]
        if X2ARRAY!=None:
                params[vardict[scanvar2]]=X2ARRAY[i]
        return DR_get_root(i,params,start,method,scan_options,XARRAY,X2ARRAY,OMEGA,GAMMA,BY,DN,shift)

#This function allows scans around the previous position to find the root's new position.
#scandirec defines an angle along which the new starting points will be aligned
def next_start(X, F, G, i, j,scandirec=0,radius=1.,method='default'):
        radius*=np.sqrt(abs(F[i-1]**2+G[i-1]**2))
        df=j*radius*np.cos(scandirec)
        dg=j*radius*np.sin(scandirec)
        if 'default' in method: #use old root position as starting point, plus potential scandirec modifications
                fval=F[i]+df
                gval=G[i]+dg
        elif 'predict' in method: #use first order prediction of the root's new position
                if i>1:
                        p=np.polyfit(X[i-2:i+1],F[i-2:i+1],2)
                        fval=p[0]*X[i+1]**2+p[1]*X[i+1]+p[2]+df
                        p=np.polyfit(X[i-2:i+1],G[i-2:i+1],2)
                        gval=p[0]*X[i+1]**2+p[1]*X[i+1]+p[2]+dg                     
                elif i==1:
                        fval=F[i]+np.diff(F)[i-1]/np.diff(X)[i-1]*np.diff(X)[i]+df
                        gval=G[i]+np.diff(G)[i-1]/np.diff(X)[i-1]*np.diff(X)[i]+dg
                else:
                        fval=F[i]
                        gval=G[i]
                if 'ldr' in method:
                        gval=-0.75*abs(G[i]+np.diff(G)[i-1]/np.diff(X)[i-1]*np.diff(X)[i])+dg
        elif 'ldr' in method: #lower damping rate (designed to prefer more weakly damped modes)
                gval=-0.75*abs(G[i])+dg #always aim for lower damping rate so the less damped mode will be chosen
                fval=F[i]+0.03*F[i]*np.sign(np.diff(X)[-1])+df #include slight frequency shift (depending on the scan direction) in order to help jump across resonances
        #lower limit for damping rates, to prevent runaway to extremely small values 
        if abs(gval)<1e-13: gval=np.sign(gval)*1e-13
        return [fval,gval] 

def DR_get_root(index=0,params=[],start=[1.0,0.0],method='default',scan_options=[1.,1000.,1.],
                XARRAY=[],X2ARRAY=[],OMEGA=[],GAMMA=[],BY=[],DN=[],shift=0):
	# finds the root by using different starting postions
	# as well as checking to make sure the root it converged to 
	# is the correct one.
        radius,dev_limit,cp_limit=scan_options

        if index==0:
                omega,gamma,By,dn,err,res=DR_Solve(params,start)
                while omega is None:
                        start[0]+=start[0]/10**5 #might not work great.
                        omega,gamma,By,dn,err,res=DR_Solve(params,start)
                return [omega,gamma,By,dn,err]

        def Scan(jstart,jstep,jmax,scandirec=0,radius=1.0,probe=False):
                j=jstart
                while (j<=jmax): 
                        start=next_start(XARRAY,OMEGA,GAMMA,index-1-shift,j,scandirec,radius,method)
                        print("start =",start,file=outfile)

                        if probe:
                                return DR_point(params,start)
                        else:
                                omega,gamma,By,dn,err,res0=DR_Solve(params,start)

                        if err==1 and (By is not None): #no error
                                if abs(omega)<1e6 and abs(gamma)<1e6: 
					#to prevent non-moving solutions and ones that go crazy.
                                        if index==0 or abs(omega-OMEGA[index-1-shift])>1e-9 or abs(gamma-GAMMA[index-1-shift])>1e-9:
                                                limit=np.float64(dev_limit)
                                                cpl=cp_limit*2*np.pi
                                                omega_comp=OMEGA[index-1-shift]
                                                gamma_comp=GAMMA[index-1-shift]
                                                d1=abs(log10(abs(omega))-log10(abs(omega_comp)))
                                                d2=0.
                                                d3=abs(log10(abs(By))-log10(abs(BY[index-1-shift])))
                                                d4=abs(log10(abs(dn))-log10(abs(DN[index-1-shift])))
                                                d5=abs(np.angle(By)-np.angle(BY[index-1-shift]))
                                                d6=abs(np.angle(dn)-np.angle(DN[index-1-shift]))
                                                if (all([d3,d4]<=limit) and (d1<=limit or abs(omega) < 1e-3) and (d2<=limit or abs(gamma) < 1e-3)) and (d5<cpl or abs(d5-2*np.pi)<cpl) and (d6<cpl or abs(d6-2*np.pi)<cpl): 
							# very close to guess
							# allow gamma, when really small, some more leeway.
                                                        return [omega,gamma,By,dn,err]
                                                else:
							#check for entropy mode
                                                        deltax=log10(XARRAY[index-1-shift])-log10(XARRAY[index-2-shift])
                                                        if abs(OMEGA[index-2-shift])!=0.:
                                                                deltay=log10(abs(OMEGA[index-1-shift]))-log10(abs(OMEGA[index-2-shift]))
                                                        else:
                                                                deltay=log10(abs(OMEGA[index-1-shift]))
                                                        slope=deltay/deltax
                                                        if (abs(omega) < 1e-3 and (slope < -2. or abs(OMEGA[index-1-shift])<1e-3) 
							    and (d2<limit or abs(gamma) < 1e-3)):
                                                                return [omega,gamma,By,dn,err]
                                                        d_arr=[d1>limit,d2>limit,d3>limit,d4>limit,not (d5<cpl or abs(d5-2*np.pi)<cpl),not (d6<cpl or abs(d6-2*np.pi)<cpl)]
                                                        if any(d_arr[0:2]): print('Discontinuity in omega/gamma, %7.4e, %7.4e' %(d1, d2),file=outfile)
                                                        if any(d_arr[2:4]): print('Discontinuity in amplitudes, %7.4e, %7.4e' %(d3, d4),file=outfile)
                                                        if any(d_arr[4:]) : print('Discontinuity in cross phases, %7.4e, %7.4e' %(d5, d6),file=outfile)
                                                        update_output('\r Continuity problem on index %d, attempting to recover mode... ' %(index+1))
                                                        j+=jstep
                                        else:
						#print("too close to start")
                                                j+=jstep
                                else:
					#print ("too far away")
                                        j+=jstep
                        else:
                                print("Eigenvalue smallness condition not fulfilled, or solver error",file=outfile)
                                j+=jstep
                return None

        # Run many scans to increase likelihood of finding root.
        scan_1=Scan(0,1,0,scandirec=0.) #original settings: 0,1,1,0!
        if scan_1 is None:
                best_res=1e10
                steps=[0.1,0.2,0.5,1.0]
                for step in steps:
                        if not scan_1==None: break
                        scandirec=-np.pi
                        while (scan_1 is None and scandirec<np.pi):
                                #We don't perform a full root solve, but only probe the dispersion relation
                                #around the starting point. We only perform a root solve for the scan direction
                                #which returned the best residual, then we proceed to the next probing radius.
                                res=Scan(-1,1,-1,scandirec,radius=step*radius,probe=True)
                                best_scandirec=scandirec
                                if res<best_res: 
                                        best_res=res
                                        best_scandirec=scandirec
                                scandirec+=np.pi/4
                        scan_1=Scan(-1,1,-1,best_scandirec,radius=step*radius)
                if scan_1 is None:
                        return None,None,None,None,0
        return scan_1

