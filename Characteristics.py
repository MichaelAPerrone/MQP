import numpy as np
import matplotlib.pyplot as plt
from matplotlib import cm
from scipy import ndimage
from matplotlib.colors import LightSource
from matplotlib import cbook

import matplotlib as mpl
from mpl_toolkits.mplot3d import Axes3D
#from mayavi import mlab

#Set the material phase velocities in materials 1 and 2
#a1=.55;a2=1.1; 
a1=1.0;a2=0.8; 

#Set the checkerboard geometrical paramaters
tau=1.0;n1=.5;
eps=1.0;m1=.5;

#Set how many spatial-temporal periods the characteristics will be plotted over
#  as well as setting the number of times 
p_T0=0;p_Tf=10;Nt=1000;
p_a =0;p_b =10;Nx=100;

#Sets the maximum amount of smoothing allowed.
Pe=min([m1*eps,(1-m1)*eps,n1*tau,(1-n1)*tau])/2

#Set the amount of smoothing to be used in FG material
alp=(.001)*Pe;bet=(.001)*Pe;

alpha=alp;beta =bet;

#Set the boundary of the checkerboard corresponding to the above
T0=p_T0*tau; Tf=p_Tf*tau;
a =p_T0*eps;  b=p_b*eps;


t_stor=np.linspace(T0,Tf,Nt);
X,T = np.mgrid[a:b:Nx*1j, T0:Tf:Nt*1j]; dx=X[1,0]-X[0,0]
Xmod,Tmod = np.mgrid[0.0:eps:Nx*1j, 0.0:tau:Nt*1j];

T0_c=0*tau; Tf_c=Tf;
#a_c=alpha+.01;b_c=m1*eps-alpha*eps-.07;

a_c=-(1-m1)*eps;b_c=eps+.3*m1*eps;
N_c=20;

X_c,T_c = np.mgrid[a_c:b_c:N_c*1j, T0_c:Tf_c:Nt*1j];
dt=T_c[0,1]-T_c[0,0];

##def LC_Loc1(delta,tau,a1,a2,m,n,alpha,beta):
##  lb=a2/a1;
##  Orig=(-lb)*delta+a2*tau+(1-lb)*lb*m*delta+(-1+lb)*a2*n*tau;
##  
##  New= beta*((-3/4)*a2+(1/4)*a1+(1/4)*a2*lb+(1/4)*a2*lb**2)
##  return ((Orig+New)/(1.0-lb**2))
def LC_Loc1(delta,tau,a1,a2,m,n,alpha,beta):
  lb=a2/a1;
  C1=(-lb**2+2*lb*np.log(lb)+1.0)/(lb**2-2*lb+1.0);
  C2=a2/4*(-1.0+1/lb)
  C3=lb/(lb+1.0)
  C4=-a2/(lb+1)
  C5=(-a2+lb)/(lb**2-1)
  
  New= C1*alpha+C2*beta+C3*m+C4*m+C5
  return New


#################################Definition of Material Geometry################
############ Standard Checkerboard
def f_u(x,t,u1,u2):
    def M1(x):
        return (u1*(np.mod(x,eps)<m1*eps) + u2*(np.mod(x,eps)>=m1*eps));
    def M2(x):
        return (u2*(np.mod(x,eps)<m1*eps) + u1*(np.mod(x,eps)>=m1*eps));
    return (M1(x)*(np.mod(t,tau)<n1*tau)+M2(x)*(np.mod(t,tau)>=n1*tau))

############ Functionallay Graded Checkerboard (tanh function)
####alpha=.75;
####beta=.75;

#alpha=0.4;
#beta=0.4;

#def f_u(x,t,u1,u2):
    #def smt(x,t):
        #return np.tanh(np.sin(2*np.pi*x)/alpha)*np.tanh(np.sin(2*np.pi*t)/beta)
    #return (u1+u2)/(2.0)+(u1-u2)/2.0*smt(x,t)
 
def l(xi,xi1,xi2,y1,y2):
    m=(y2-y1)/(xi2-xi1);
    return y1+m*(xi-xi1)
    
def p(xi,eta,xi1,xi2,eta1,eta2,y1,y2):
    mz=(y2 - y1 )/(xi2 -xi1 );
    mt=(y2 - y1 )/(eta2-eta1);
    return y1+mz*(xi-xi1)+mt*(eta-eta1);

#############One ST inclusion

#def f_u(x,t,u1,u2):
    #Wx=eps/10;
    #Wt=tau/5;
    #xTi=np.mod(x,eps);
    #tTi=np.mod(t,tau);
    #z1=m1*eps-Wx;
    #z2=m1*eps+Wx;
    #t1=n1*tau-Wt;
    #t2=n1*tau+Wt;
    
    ##def M1(x):
        ##return u1*(xTi < z1) + l(x,z1,z2,u1,u2)*((z1<=xTi)*(xTi < z2))+u2*(z2<= xTi)
    ##def M2(x,t):
        ##return (l(t,t1,t2,u1,u2)*(xTi < z1) + 
	       ##p(x,t,z1,z2,t1,t2,u1,u2)*((z1 <= xTi)*(xTi < z2)*(xTi-z1 <((z2-z1)/(t1-t2))*(tTi-t2)))+
	       ##p(x,t,z2,z1,t2,t1,u1,u2)*((z1 <= xTi)*(xTi < z2)*(xTi-z1>=((z2-z1)/(t1-t2))*(tTi-t2)))+
	       ##l(t,t1,t2,u2,u1)*(z2 <= xTi))
    ##def M3(x):
        ##return u2*(xTi < z1) + l(x,z2,z1,u1,u2)*((z1<=xTi)*(xTi < z2))+u1*(z2<= xTi);
    
    #def M1(x):
        #return u1*(x < z1) + l(x,z1,z2,u1,u2)*((z1<=x)*(x < z2))+u2*(z2<= x)
    #def M2(x,t):
        #return (l(t,t1,t2,u1,u2)*(x < z1) + 
	       #p(x,t,z1,z2,t1,t2,u1,u2)*((z1 <= x)*(x < z2)*(x-z1 <((z2-z1)/(t1-t2))*(t-t2)))+
	       #p(x,t,z2,z1,t2,t1,u1,u2)*((z1 <= x)*(x < z2)*(x-z1>=((z2-z1)/(t1-t2))*(t-t2)))+
	       #l(t,t1,t2,u2,u1)*(z2 <= x))
    #def M3(x):
        #return u2*(x < z1) + l(x,z2,z1,u1,u2)*((z1<=x)*(x < z2))+u1*(z2<= x);
    
    #def Fin(x,t):
	#return M1(x)*(t <t1)+M2(x,t)*((t1 <= t)*(t <t2))+M3(x)*(t2<= t)
    
    #return Fin(xTi,tTi)


###########################Linear FG Checkerboard    
#def py(x,z1,z2,t,t1,t2,u1,u2):
#    return (p(x,t,z1,z2,t1,t2,u1,u2)*((z1 <= x)*(x < z2)*(x-z1 <((z2-z1)/(t1-t2))*(t-t2)))+
#	    p(x,t,z2,z1,t2,t1,u1,u2)*((z1 <= x)*(x < z2)*(x-z1>=((z2-z1)/(t1-t2))*(t-t2))))
    
#Wx=alpha
#Wt=beta;    
#def f_u(x,t,u1,u2):
#    xTi=np.mod(x,eps);
#    tTi=np.mod(t,tau);
#commented out so that non functionally graded checkerboard can be used - to use functionally graded code need to uncomment this and comment out f_u above. Namespaces man...    
    #def M1(x):
        #return u1*(xTi < z1) + l(x,z1,z2,u1,u2)*((z1<=xTi)*(xTi < z2))+u2*(z2<= xTi)
    #def M2(x,t):
        #return (l(t,t1,t2,u1,u2)*(xTi < z1) + 
	       #p(x,t,z1,z2,t1,t2,u1,u2)*((z1 <= xTi)*(xTi < z2)*(xTi-z1 <((z2-z1)/(t1-t2))*(tTi-t2)))+
	       #p(x,t,z2,z1,t2,t1,u1,u2)*((z1 <= xTi)*(xTi < z2)*(xTi-z1>=((z2-z1)/(t1-t2))*(tTi-t2)))+
	       #l(t,t1,t2,u2,u1)*(z2 <= xTi))
    #def M3(x):
        #return u2*(xTi < z1) + l(x,z2,z1,u1,u2)*((z1<=xTi)*(xTi < z2))+u1*(z2<= xTi);
    
#    def M0(x,t):      
#        return (py(x,-Wx,Wx,t,-Wt,Wt,u1,u2)*(x<Wx)+
#                l(t,-Wt,Wt,u2,u1)*(x>=Wx)*(x<m1*eps-Wx)+
#                py(x,m1*eps-Wx,m1*eps+Wx,t,-Wt,Wt,u2,u1)*(x>=m1*eps-Wx)*(x<m1*eps+Wx)+ 
#		l(t,-Wt,Wt,u1,u2)*(m1*eps+Wx<=x)*(x<eps-Wx)+
#                py(x,eps-Wx,eps+Wx,t,-Wt,Wt,u1,u2)*(x>eps-Wx))
    
#    def M1(x,t,u1,u2):
#        z1=m1*eps-Wx;
#        z2=m1*eps+Wx;
#        return (l(x,-Wx,Wx,u2,u1)*(x<Wx)+
#                u1*(Wx<=x)*(x < z1) + 
#                l(x,z1,z2,u1,u2)*((z1<=x)*(x < z2))+
#                u2*(z2<= x)*(x<eps-Wx)+
#                l(x,eps-Wx,eps+Wx,u2,u1)*(eps-Wx<=x))

#    def M2(x,t):
#        t1=n1*tau-Wt;
#        t2=n1*tau+Wt;
#        z1=m1*eps-Wx;
#        z2=m1*eps+Wx;
#        return (py(x,-Wx,Wx,t,t1,t2,u2,u1)*(x<Wx)+
#	       l(t,t1,t2,u1,u2)*(Wx<=x)*(x < z1) + 
#	       p(x,t,z1,z2,t1,t2,u1,u2)*((z1 <= x)*(x < z2)*(x-z1 <((z2-z1)/(t1-t2))*(t-t2)))+
#	       p(x,t,z2,z1,t2,t1,u1,u2)*((z1 <= x)*(x < z2)*(x-z1>=((z2-z1)/(t1-t2))*(t-t2)))+
#	       l(t,t1,t2,u2,u1)*(z2 <= x)*(x<eps-Wx)+
#	       py(x,eps-Wx,eps+Wx,t,n1*tau-Wt,n1*tau+Wt,u2,u1)*(eps-Wx<=x))
#    def M3(x,t):
#        z1=m1*eps-Wx;
#        z2=m1*eps+Wx;
#        return l(x,-Wx,Wx,u1,u2)*(x<Wx)+u2*(Wx<=x)*(x < z1) + l(x,z2,z1,u1,u2)*((z1<=x)*(x < z2))+u1*(z2<= x)*(x<eps-Wx)+l(x,eps-Wx,eps+Wx,u1,u2)*(eps-Wx<=x);
#    def M3(x,t):
#        return M1(x,t,u2,u1)
    
#    def M4(x,t):
#        return (py(x,-Wx,Wx,t,tau-Wt,tau+Wt,u1,u2)*(x<Wx)+
#                 l(t,tau-Wt,tau+Wt,u2,u1)*(x>=Wx)*(x<m1*eps-Wx)+
#                py(x,m1*eps-Wx,m1*eps+Wx,t,tau-Wt,tau+Wt,u2,u1)*(x>=m1*eps-Wx)*(x<m1*eps+Wx)+ 
#		l(t,tau-Wt,tau+Wt,u1,u2)*(m1*eps+Wx<=x)*(x<eps-Wx)+py(x,eps-Wx,eps+Wx,t,tau-Wt,tau+Wt,u1,u2)*(x>=eps-Wx))
      
#    def Fin(x,t):
#        t1=n1*tau-Wt;
#        t2=n1*tau+Wt;
#	return (M0(x,t)*(t<Wt)+
#                M1(x,t,u1,u2)*(Wt<=t)*(t <t1)+
#                M2(x,t)*((t1 <= t)*(t <t2))+
#                M3(x,t)*(t2<= t)*(t<tau-Wt)+
#                M4(x,t)*(tau-Wt<=t))

#    return Fin(xTi,tTi)

#def f_u(x,t,u1,u2):
    #Wx=eps/10;
    #Wt=tau/10;
    #xTi=np.mod(x,eps);
    #tTi=np.mod(t,tau);
    #z1=m1*eps-Wx;
    #z2=m1*eps+Wx;
    #t1=n1*tau-Wt;
    #t2=n1*tau+Wt;
    
    #def pyr(x,y,x1,x2,y1,y2,u1,u2):
        #m=(y2-y1)/(x1-x)
        #return (p(x,y,x1,x2,y1,y2,u1,u2)*((x-x1)< m*(y-y1))+
		#p(x,y,x2,x1,y2,y1,u1,u2)*((x-x1)>=m*(y-y1)))
       
    
    
    #def M1(x):
        #return l(x,Wx,-Wx,u1,u2)*(x<Wx)+u1*((Wx<=x)*(x < z1)) + l(x,z1,z2,u1,u2)*((z1<=x)*(x < z2))+u2*((z2<= x)*(x<eps-Wx))+l(x,eps-Wx,eps+Wx,u2,u1)*(eps-Wx<=x)
    #def M2(x,t):
        #return (p(x,t,-Wx,Wx,t1,t2,u2,u1)*((x<Wx)*(x+Wx< ((2*Wx)/(t1-t2))*(t-t2))) +
		#p(x,t,Wx,-Wx,t2,t1,u2,u1)*((x<Wx)*(x+Wx>=((2*Wx)/(t1-t2))*(t-t2)))+
	        #l(t,t1,t2,u1,u2)*((Wx<=x)*(x < z1)) + 
	        #p(x,t,z1,z2,t1,t2,u1,u2)*((z1 <= x)*(x < z2)*(x-z1 <((z2-z1)/(t1-t2))*(t-t2)))+
	        #p(x,t,z2,z1,t2,t1,u1,u2)*((z1 <= x)*(x < z2)*(x-z1>=((z2-z1)/(t1-t2))*(t-t2)))+
	        #l(t,t1,t2,u2,u1)*((z2<=x)*(x < eps-Wx)) + 
	        #pyr(x,t,eps-Wx,eps+Wx,t1,t2,u1,u2)*(eps-Wx<=x))
    #def M3(x):
        #return l(x,Wx,-Wx,u2,u1)*(x<Wx)+u2*((Wx<x)*(x < z1)) + l(x,z2,z1,u1,u2)*((z1<=x)*(x < z2))+u1*((z2<= x)*(x<eps-Wx))+l(x,eps-Wx,eps+Wx,u1,u2)*(eps-Wx<=x);
    
    #def Fin(x,t):
	#return M1(x)*(t <t1)+M2(x,t)*((t1 <= t)*(t <t2))+M3(x)*(t2<= t)
    
    #return Fin(xTi,tTi)

###################### Solution of Characteristic Equation
C_c=np.zeros(X_c.shape);

#### Rk4 Method

N_unstable=2;
X_unstable=np.zeros((N_unstable,X_c.shape[-1]));
T_unstable=np.zeros(X_unstable.shape);

N_stable=2;
X_stable=np.zeros((N_stable,X_c.shape[-1]));
T_stable=np.zeros(X_stable.shape);

T_unstable[:,0]=T0+T_unstable[:,0];
Unstable_init=m1*eps+LC_Loc1(eps,tau,a2,a1,1-m1,n1,alpha,beta);
X_unstable[:,0]=np.array([Unstable_init,Unstable_init+eps]);

T_stable[:,0]=T0+T_stable[:,0];
Stable_init=LC_Loc1(eps,tau,a1,a2,m1,n1,alpha,beta);
X_stable[:,0]=np.array([Stable_init,Stable_init+eps]);

C_avg=np.zeros(X_c[:,0].shape)

Np=np.int((5.0/6)*Nt);

for n in range(Nt-1):
#C_c[:,n]=f_u(X_c[:,n],T_c[:,n],a1,a2);
    k1=f_u(X_c[:,n]        ,T_c[:,n]     ,a1,a2);						
    k2=f_u(X_c[:,n]+k1*dt/2,T_c[:,n]+dt/2,a1,a2);	
    k3=f_u(X_c[:,n]+k2*dt/2,T_c[:,n]+dt/2,a1,a2);	
    k4=f_u(X_c[:,n]+k3*dt  ,T_c[:,n]+dt  ,a1,a2);
    X_c[:,n+1]=X_c[:,n]+(dt/6)*(k1+2*k2+2*k3+k4);
	
    k1_u=f_u(X_unstable[:,n]          ,T_unstable[:,n]     ,a1,a2); 
    k2_u=f_u(X_unstable[:,n]+k1_u*dt/2,T_unstable[:,n]+dt/2,a1,a2);
    k3_u=f_u(X_unstable[:,n]+k2_u*dt/2,T_unstable[:,n]+dt/2,a1,a2);
    k4_u=f_u(X_unstable[:,n]+k3_u*dt  ,T_unstable[:,n]+dt  ,a1,a2);
    X_unstable[:,n+1]=X_unstable[:,n]+(dt/6)*(k1_u+2*k2_u+2*k3_u+k4_u);T_unstable[:,n+1]=T_unstable[:,n]+dt;

    k1_s=f_u(X_stable[:,n]          ,T_stable[:,n]     ,a1,a2); 
    k2_s=f_u(X_stable[:,n]+k1_s*dt/2,T_stable[:,n]+dt/2,a1,a2);
    k3_s=f_u(X_stable[:,n]+k2_s*dt/2,T_stable[:,n]+dt/2,a1,a2);
    k4_s=f_u(X_stable[:,n]+k3_s*dt  ,T_stable[:,n]+dt  ,a1,a2);
    X_stable[:,n+1]=X_stable[:,n]+(dt/6)*(k1_s+2*k2_s+2*k3_s+k4_s);T_stable[:,n+1]=T_stable[:,n]+dt;
	
    if t_stor[n]>=t_stor[Np]:
	C_avg=C_avg+f_u(X_c[:,n],T_c[:,n],a1,a2)
	
C_avg=C_avg/np.size(t_stor[Np:-1])	
	
#### diff1=np.abs(X_c[:, 0][1]-X_c[:, 0][0])
#### diff2=np.abs(X_c[:,-1][1]-X_c[:,-1][0])

#### TotalDiff=diff2-diff1;

Np=np.int((2.0/3)*C_c.shape[1]);

StorAvg=np.zeros(C_c.shape[0]);

#for n in range(C_c.shape[0]):
#    StorAvg[n]=np.average(C_c[4,Np:-1])

CharNumber=0;
Xcmod=np.mod(X_c,eps);Tcmod=np.mod(T_c,tau);

## Create Custom ColorMap Function
def make_cmap(colors, position=None, bit=False):
    '''
    make_cmap takes a list of tuples which contain RGB values. The RGB
    values may either be in 8-bit [0 to 255] (in which bit must be set to
    True when called) or arithmetic [0 to 1] (default). make_cmap returns
    a cmap with equally spaced colors.
    Arrange your tuples so that the first color is the lowest value for the
    colorbar and the last is the highest.
    position contains values from 0 to 1 to dictate the location of each color.
    '''
    import matplotlib as mpl
    import numpy as np
    bit_rgb = np.linspace(0,1,256)
    if position == None:
        position = np.linspace(0,1,len(colors))
    else:
        if len(position) != len(colors):
            sys.exit("position length must be the same as colors")
        elif position[0] != 0 or position[-1] != 1:
            sys.exit("position must start with 0 and end with 1")
    if bit:
        for i in range(len(colors)):
            colors[i] = (bit_rgb[colors[i][0]],
                         bit_rgb[colors[i][1]],
                         bit_rgb[colors[i][2]])
    cdict = {'red':[], 'green':[], 'blue':[]}
    for pos, color in zip(position, colors):
        cdict['red'].append((pos, color[0], color[0]))
        cdict['green'].append((pos, color[1], color[1]))
        cdict['blue'].append((pos, color[2], color[2]))

    cmap = mpl.colors.LinearSegmentedColormap('my_colormap',cdict,256)
    return cmap
  
RedCharCmap=make_cmap([(220,220,220),(140,0,26)],bit=True)

colors = [(100,200,255),(255,255,150)]
MyCmap=make_cmap(colors,bit=True)
  
########## Figure-1 ################# 

plt.figure(1)
gr=(0.863,0.863,0.863)

plt.pcolormesh(X,T,f_u(X,T,a1,a2),cmap=MyCmap)
plt.colorbar()

n_Contours=40;

#plt.contour(X,T,f_u(X,T,a1,a2),contours=n_Contours,colors='k')


plt.xlim([0.0,1.0*eps])
#plt.xlim([0.0,1.0*eps+m1*eps])
plt.ylim([0.0,1.0*tau])


plt.title(r'$\delta =$'+str(eps) +  \
	  r', $\tau =$' + str(tau) + ', $m$ =' + str(m1) + ', $n$=' +str(n1) +\
	  r', $\alpha_1=$'   + str(a1)  + r', $\alpha_2$=' + str(a2)+", \n"+
	   r'$p$='+str(alpha)+r', $q=$'+str(beta)+r',    $dt=$'+str(dt))

plt.ylabel('$t\in$['+str(p_T0)+r'$\tau, $'    +str(p_Tf)+r'$\tau$]')
plt.xlabel('$z\in$['+str(p_a )+r'$\epsilon, $'+str(p_b) +r'$\epsilon$]')

frame1=plt.gca()

frame1.axes.get_xaxis().set_ticks([])
frame1.axes.get_yaxis().set_ticks([])

plt.plot(np.rot90(X_c),np.rot90(T_c),'k-',linewidth=2.0,color=gr);

plt.savefig("Checkerboard.png",dpi=1000)

plt.clf()
########## Figure-2 ################# 

plt.figure(2)
gr=(0.863,0.863,0.863)

plt.pcolormesh(X,T,f_u(X,T,a1,a2),cmap=MyCmap)
plt.colorbar()

n_Contours=40;

#plt.contour(X,T,f_u(X,T,a1,a2),contours=n_Contours,colors='k')

#This is a hack: 0.5 should be m and n, but for the case I am working with m=n=0.5 so it happens to work
plt.xlim([0.5*eps,1.5*eps])
#plt.xlim([0.0,1.0*eps+m1*eps])
plt.ylim([0.5*tau,1.5*tau])


plt.title(r'$\delta =$'+str(eps) +  \
	  r', $\tau =$' + str(tau) + ', $m$ =' + str(m1) + ', $n$=' +str(n1) +\
	  r', $\alpha_1=$'   + str(a1)  + r', $\alpha_2$=' + str(a2)+", \n"+
	   r'$p$='+str(alpha)+r', $q=$'+str(beta)+r',    $dt=$'+str(dt))

plt.ylabel('$t\in$['+str(p_T0)+r'$\tau, $'    +str(p_Tf)+r'$\tau$]')
plt.xlabel('$z\in$['+str(p_a )+r'$\epsilon, $'+str(p_b) +r'$\epsilon$]')

frame1=plt.gca()

frame1.axes.get_xaxis().set_ticks([])
frame1.axes.get_yaxis().set_ticks([])

plt.plot(np.rot90(X_c),np.rot90(T_c),'k-',linewidth=2.0,color=gr);

plt.savefig("Checkerboard2.png",dpi=1000)

plt.clf()
#Added to plot second family

plt.figure(3)

plt.pcolormesh(X,T,f_u(X,T,a1,a2),cmap=MyCmap)
plt.colorbar()

n_Contours=1;

#plt.contour(X,T,f_u(X,T,a1,a2),colors='k',linestyle='dotted')

plt.xlim([0.0,b])
plt.ylim([0.0,Tf])


plt.title(r'$\delta =$'+str(eps) +  \
	  r', $\tau =$' + str(tau) + r', $m$ =' + str(m1) + r', $n$=' +str(n1) +\
	  r', $\alpha_1=$'   + str(a1)  + r', $\alpha_2$=' + str(a2)+", \n"+
	   r'$p$='+str(alpha)+r', $q=$'+str(beta)+r',    $dt=$'+str(dt))

plt.ylabel('$t\in$['+str(p_T0)+r'$\tau, $'    +str(p_Tf)+r'$\tau$]')
plt.xlabel('$z\in$['+str(p_a )+r'$\epsilon, $'+str(p_b) +r'$\epsilon$]')

frame1=plt.gca()

frame1.axes.get_xaxis().set_ticks([])
frame1.axes.get_yaxis().set_ticks([])

#Char=plt.plot(np.rot90(X_c),np.rot90(T_c),colormap=RedCharCmap,linewidth=1.0);

p_R=Nt/p_Tf;n_R=p_Tf-2;

gr=(0.863,0.863,0.863)
re=(0.941,0.0,0.102)

Char=plt.plot(np.rot90(X_c[:,0:n_R*p_R]),np.rot90(T_c[:,0:n_R*p_R]),color=gr,linewidth=1.0);

Char=plt.plot(np.rot90(X_c[:,n_R*p_R:Nt-1]),np.rot90(T_c[:,n_R*p_R:Nt-1]),color=re,linewidth=1.0);
#plt.scatter(np.rot90(X_c),np.rot90(T_c))

plt.savefig("Checkerboard_Full.png",dpi=1000)

plt.clf()



plt.figure(4)

colors = [(100,200,255),(255,255,150)]
MyCmap=make_cmap(colors,bit=True)

plt.pcolormesh(X,T,f_u(X,T,a1,a2),cmap=MyCmap)
plt.colorbar()

n_Contours=2;

#plt.contour(X,T,f_u(X,T,a1,a2),contours=n_Contours,colors='k')

plt.xlim([9*eps,10*eps])
plt.ylim([9*tau,10*tau])

plt.title(r'$\delta =$'+str(eps) +  \
	  r', $\tau =$' + str(tau) + ', $m$ =' + str(m1) + ', $n$=' +str(n1) +\
	  r', $\alpha_1=$'   + str(a1)  + r', $\alpha_2$=' + str(a2)+", \n"+
	   r'$p$='+str(alpha)+r', $q=$'+str(beta)+r',    $dt=$'+str(dt))

plt.ylabel('$t\in$['+str(4)+r'$\tau, $'    +str(5)+r'$\tau$]')
plt.xlabel('$z\in$['+str(4 )+r'$\delta, $'+str(5) +r'$\epsilon$]')

frame1=plt.gca()

frame1.axes.get_xaxis().set_ticks([])
frame1.axes.get_yaxis().set_ticks([])

plt.plot(np.rot90(X_c),np.rot90(T_c),'-',color=re,linewidth=2.0);

plt.savefig("Checkerboard_LC.png",dpi=1000)

plt.clf()





plt.figure(4)

colors = [(100,200,255),(255,255,150)]
MyCmap=make_cmap(colors,bit=True)

plt.pcolormesh(X,T,f_u(X,T,a1,a2),cmap=MyCmap)
plt.colorbar()

n_Contours=2;

#plt.contour(X,T,f_u(X,T,a1,a2),contours=n_Contours,colors='k')

plt.xlim([0.0*eps,3.0*eps+m1*eps])
plt.ylim([0.0*tau,3.0*tau])

plt.title('Checkerboard for '+r'$\delta =$'+str(eps) +  \
	  r', $\tau =$' + str(tau) + ', $m$ =' + str(m1) + ', $n$=' +str(n1) +\
	  r', $\alpha_1=$'   + str(a1)  + r', $\alpha_2$=' + str(a2)+", \n"+
	   r'$p$='+str(alpha)+r', $q=$'+str(beta))

plt.ylabel('$t\in$['+str(0)+r'$\tau, $'    +str(3)+r'$\tau$]')
plt.xlabel('$z\in$['+str(0 )+r'$\delta, $'+str(3) +r'$\delta$]')

frame1=plt.gca()

frame1.axes.get_xaxis().set_ticks([])
frame1.axes.get_yaxis().set_ticks([])

#plt.plot(np.rot90(X_c),np.rot90(T_c),'-',color='darkred',linewidth=2.0);

plt.savefig("Checkerboard_NoChar.png",dpi=1000)

plt.clf()

############# Figure-2 ################# 

###X_unstable_mod=np.mod(X_unstable,eps);T_unstable_mod=np.mod(T_unstable,tau);

###plt.figure(2)


###plt.pcolormesh(Xmod,Tmod,f_u(Xmod,Tmod,a1,a2),cmap=MyCmap)
###plt.colorbar()

###n_Contours=20;

###plt.contour(Xmod,Tmod,f_u(Xmod,Tmod,a1,a2),contours=n_Contours)


####plt.xlim([0.0,eps])
####plt.ylim([0.0,tau])
####plt.colorbar()

####plt.title('Characteristics for '+r'$\epsilon =$'+str(eps) +  \
	  ####r', $\tau =$' + str(tau) + ', $m$ =' + str(m1) + ', $n$=' +str(n1) +\
	  ####', $a_1=$'   + str(a1)  + ', $a_2$=' + str(a2)+", \n"+
	   ####r'$\alpha$='+str(alpha)+r', $\beta=$'+str(beta))

####plt.ylabel('$t\in$['+str(p_T0)+r'$\tau, $'    +str(p_Tf)+r'$\tau$]')
####plt.xlabel('$z\in$['+str(p_a )+r'$\epsilon, $'+str(p_b) +r'$\epsilon$]')

###frame1=plt.gca()

###frame1.axes.get_xaxis().set_ticks([])
###frame1.axes.get_yaxis().set_ticks([])


####frame1.axes.get_xaxis().set_visible(False)
####frame1.axes.get_yaxis().set_visible(False)

#######XcRot=np.rot90(Xcmod[CharNumber,:]);
#######TcRot=np.rot90(Tcmod[CharNumber,:]);

###plt.plot(Xcmod[CharNumber,:],Tcmod[CharNumber,:],'b*',linewidth=1.0);
###plt.plot(X_unstable_mod[0,:],T_unstable_mod[0,:],'g*',linewidth=3.0);



####plt.plot(XcRot,TcRot,'b*',linewidth=1.0);

##### plt.plot(np.rot90(X_unstable),np.rot90(T_unstable),'k--',linewidth=3.0);
####plt.plot(np.rot90(X_stable),np.rot90(T_stable),'k-',linewidth=3.0);

####plt.savefig("Checkerboard.png",dpi=1000)

####plt.savefig('2Characteristics for '+r'$\epsilon =$'+str(eps) +  \
	  ####r', $\tau =$' + str(tau) + ', $m$ =' + str(m1) + ', $n$=' +str(n1) +\
	  ####', $a_1=$'   + str(a1)  + ', $a_2$=' + str(a2) \
	    ####+r"$\alpha$"+ str(alpha)+r"$\beta$"+str(beta) \
	    ####+", p_T0="+str(p_T0)+"P_Tf"+str(p_Tf)  \
	    ####+".png",dpi=100)

#######plt.savefig('3Characteristics for '+r'$\epsilon =$'+str(eps) +  \
	  #######r', $\tau =$' + str(tau) + ', $m$ =' + str(m1) + ', $n$=' +str(n1) +\
	  #######', $a_1=$'   + str(a1)  + ', $a_2$=' + str(a2) \
	    #######+r"$\alpha$"+ str(alpha)+r"$\beta$"+str(beta) \
	    #######+", p_T0="+str(p_T0)+"P_Tf"+str(p_Tf)  \
	    #######+".png",dpi=100)

############# Figure-4 ################# 

###plt.figure(4)


###X_Int=np.mod(X_c,eps);T_Int=np.mod(T_c,tau);

###Tol=.01;

###X_0=np.where(X_Int == 0.0);
###T_0=np.where(T_Int == 0.0);

###Xabs1=np.where(np.abs(X_Int-eps)<Tol);
###Xabs2=np.where(np.abs(X_Int-m1*eps)<Tol);
###Tabs1=np.where(np.abs(T_Int-tau)<Tol);
###Tabs2=np.where(np.abs(T_Int-n1*tau)<Tol);



###plt.pcolormesh(X,T,f_u(X,T,a1,a2),cmap=MyCmap)
###plt.colorbar()

###n_Contours=20;

###plt.contour(X,T,f_u(X,T,a1,a2),contours=n_Contours)


###plt.xlim([a,b])
###plt.ylim([T0,Tf])
####plt.colorbar()

###plt.title('Characteristics for '+r'$\epsilon =$'+str(eps) +  \
	  ###r', $\tau =$' + str(tau) + ', $m$ =' + str(m1) + ', $n$=' +str(n1) +\
	  ###', $a_1=$'   + str(a1)  + ', $a_2$=' + str(a2)+", \n"+
	   ###r'$\alpha$='+str(alpha)+r', $\beta=$'+str(beta))

####plt.ylabel('$t\in$['+str(p_T0)+r'$\tau, $'    +str(p_Tf)+r'$\tau$]')
####plt.xlabel('$z\in$['+str(p_a )+r'$\epsilon, $'+str(p_b) +r'$\epsilon$]')

###frame1=plt.gca()

###frame1.axes.get_xaxis().set_ticks([])
###frame1.axes.get_yaxis().set_ticks([])


####frame1.axes.get_xaxis().set_visible(False)
####frame1.axes.get_yaxis().set_visible(False)

#######XcRot=np.rot90(Xcmod[CharNumber,:]);
#######TcRot=np.rot90(Tcmod[CharNumber,:]);


##################Characteristics
####plt.plot(Xcmod[CharNumber,:],Tcmod[CharNumber,:],'b*',linewidth=1.0);

###plt.plot(np.rot90(X_c),np.rot90(T_c),'b-',linewidth=1.0);

###plt.plot(X_c[X_0],T_c[X_0],'r*')
###plt.plot(X_c[T_0],T_c[T_0],'r*')

###plt.plot(X_c[Xabs1],T_c[Xabs1],'k*')
###plt.plot(X_c[Xabs2],T_c[Xabs2],'k*')
###plt.plot(X_c[Tabs1],T_c[Tabs1],'g*')
###plt.plot(X_c[Tabs2],T_c[Tabs2],'g*')

####plt.plot(XcRot,TcRot,'b*',linewidth=1.0);

##### plt.plot(np.rot90(X_unstable),np.rot90(T_unstable),'k--',linewidth=3.0);
####plt.plot(np.rot90(X_stable),np.rot90(T_stable),'k-',linewidth=3.0);

####plt.savefig("Checkerboard.png",dpi=1000)

####plt.savefig('4Characteristics for '+r'$\epsilon =$'+str(eps) +  \
	  ####r', $\tau =$' + str(tau) + ', $m$ =' + str(m1) + ', $n$=' +str(n1) +\
	  ####', $a_1=$'   + str(a1)  + ', $a_2$=' + str(a2) \
	    ####+r"$\alpha$"+ str(alpha)+r"$\beta$"+str(beta) \
	    ####+", p_T0="+str(p_T0)+"P_Tf"+str(p_Tf)  \
	    ####+".png",dpi=100)


####plt.show()



####Torus1


#####R=1.0;r=0.5;


#####NN=100;
#####Theta,Phi=np.mgrid[0.0:eps:NN*1j, 0.0:tau:NN*1j];


#####X_Torus =(R+r*np.cos(2*np.pi*Phi/tau))*(np.cos(2*np.pi*Theta/eps))
#####Y_Torus =(R+r*np.cos(2*np.pi*Phi/tau))*(np.sin(2*np.pi*Theta/eps))
#####Z_Torus =   r*np.sin(2*np.pi*Phi/tau)


#####Xc_Torus =(R+r*np.cos(2*np.pi*Tcmod/tau))*(np.cos(2*np.pi*Xcmod/eps))
#####Yc_Torus =(R+r*np.cos(2*np.pi*Tcmod/tau))*(np.sin(2*np.pi*Xcmod/eps))
#####Zc_Torus =   r*np.sin(2*np.pi*Tcmod/tau)



#####mpl.rcParams['legend.fontsize'] = 10

#####fig3 = plt.figure()
#####ax3 = fig3.gca(projection='3d')


#####ax3.plot(Xc_Torus[CharNumber,:], Yc_Torus[CharNumber,:], Zc_Torus[CharNumber,:],alpha=.4)



#####G=f_u(Theta,Phi,a1,a2)

#####N=(G-G.min())/(G.max()-G.min())



#####ax3.plot_surface(X_Torus, Y_Torus, Z_Torus, cmap=MyCmap, rstride=1, cstride=1, facecolors=cm.Pastel1(N),
                       #####linewidth=0, antialiased=False,alpha=.6)


