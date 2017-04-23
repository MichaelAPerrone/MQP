#!/usr/bin/env python
# encoding: utf-8
from __future__ import division
r"""
1D polarised EM wave equation w/ variable coeficcients in space and time
==============================================

Solve variable-coefficient maxwells equations without free charges or currents, but with variable material properties in 1+1D:

.. math:: 
    \phi_t +  (\psi_x) / \epsilon(x,t) & = 0 \\ 
    \psi_t + (\phi_x) / \mu(x,t) & = 0 \\

This code utilizes a finite volume method available in Clawpack for use with the acoustics equation to simulate electromagnetic wave
propagation for a polarized plane wave in a 1+1D geometry with material properties that are variable in space and time.
Potentials are used to represent the electromagnetic system, with :math:`\phi_x & = B` and :math:`\psi_t & = B / \mu(x,t)` for magnetic field,
:math:`\psi_x & = D` and :math:`\phi_t & = D / \epsilon(x,t)` for electric displacement,
:math:`\mu(x,t) & = 1 / K` is the permeability of free space, or inverse of bulk modulus in terms of the acoustic wave equation,
and :math:`\epsilon(x,t) & = \rho(x,t)` is the permittivity of the material, or density in the acoustic wave equation.
Pressure P and Velocity V are related to the electromagnetic potentials as follows :math:`\psi & = P` :math:`\phi & = V`
This code was adapted from a 2D example, but uses only one spatial dimension
"""

# state.q 0 is equivalent to psi
# state.q 1 is equivalent to phi
# state.aux 0 is equivalent to epsilon, the permittivity
# state.aux 1 is equivalent to the wave speed. Note that in the acoustics equation rho*C = K, the bulk modulus, which is the inverse of mu
     
import numpy as np
import matplotlib.pyplot as plt

#here we define the size ax, bx, and resolution n_x of the simulated space
#ax=0.0;bx=1.0;n_x=200
ax=0.0;bx=1.0;n_x=10000 #higher resolution allows us to avoid non-physical resolution limit of power amplification for longer. Mesh refinement would be better.
SpaceStepSize = (bx - ax)/n_x #note that this is only accurate for an evenly-spaced grid. For mesh refinement an alternative must be found
ay=0.0;by=1.0;n_y=1
#n_y is just 1 because we have a 1D system, not a 2D system, though the original code was for 2D


MaterialParams = 1 #set this variable to pick which material properties we have
#Need to experiment further here with mu, epsilon, M and N, wave speed and wave impedance to get a better sense for the space and help us improve our lower bounding when considering convergence rate.

#The original code this is derived from was for sound waves, so these variables are in terms of rho and gamma, but these translate to mu and epsilon fairly easily for the equivalent 1+1D electromagnetic wave equations.
if MaterialParams == 1:
    #print('default with slight reflection\n') #This is known to work with the code that exists
    gamma=1.0;
    gamma_1=gamma;
    gamma_2=gamma+.1*gamma;
    c_1 = .6    # Sound speed (left)
    c_2 = 1.1  # Sound speed (right)
    bulk_1A = gamma_1*c_1  # Bulk modulus in left half
    bulk_2A= gamma_2*c_2  # Bulk modulus in right half
    rho_1A  = gamma_1/c_1  # Density in left half
    rho_2A  = gamma_2/c_2 # Density in right half
    bulk_1=bulk_1A
    bulk_2=bulk_2A
    rho_1=rho_1A
    rho_2=rho_2A
elif MaterialParams == 2:
    #print('No impedance mismatch or reflection\n')
    gamma=1.0;
    gamma_1=gamma;
    gamma_2=gamma;
    c_1 = .6    # Sound speed (left)
    c_2 = 1.1  # Sound speed (right)
    bulk_1A = gamma_1*c_1  # Bulk modulus in left half
    bulk_2A= gamma_2*c_2  # Bulk modulus in right half
    rho_1A  = gamma_1/c_1  # Density in left half
    rho_2A  = gamma_2/c_2 # Density in right half
    bulk_1=bulk_1A
    bulk_2=bulk_2A
    rho_1=rho_1A
    rho_2=rho_2A
elif MaterialParams == 3:
    #print('Large reflection\n')
    gamma=1.0;
    gamma_1=gamma;
    gamma_2=gamma + gamma;
    c_1 = .6    # Sound speed (left)
    c_2 = 1.1  # Sound speed (right)
    bulk_1A = gamma_1*c_1  # Bulk modulus in left half
    bulk_2A= gamma_2*c_2  # Bulk modulus in right half
    rho_1A  = gamma_1/c_1  # Density in left half
    rho_2A  = gamma_2/c_2 # Density in right half
    bulk_1=bulk_1A
    bulk_2=bulk_2A
    rho_1=rho_1A
    rho_2=rho_2A
elif MaterialParams == 4:
    #print('Large wave speed ratio\n')
    gamma=1.0;
    gamma_1=gamma;
    gamma_2=gamma;
    c_1 = .2    # Sound speed (left)
    c_2 = 1.6  # Sound speed (right)
    bulk_1A = gamma_1*c_1  # Bulk modulus in left half
    bulk_2A= gamma_2*c_2  # Bulk modulus in right half
    rho_1A  = gamma_1/c_1  # Density in left half
    rho_2A  = gamma_2/c_2 # Density in right half
    bulk_1=bulk_1A
    bulk_2=bulk_2A
    rho_1=rho_1A
    rho_2=rho_2A
elif MaterialParams == 5:
    #print('sanity check: no reflection no wavespeed so no energy change\n')
    gamma=1.0;
    gamma_1=gamma;
    gamma_2=gamma;
    c_1 = 1.0    # Sound speed (left)
    c_2 = 1.0  # Sound speed (right)
    bulk_1A = gamma_1*c_1  # Bulk modulus in left half
    bulk_2A= gamma_2*c_2  # Bulk modulus in right half
    rho_1A  = gamma_1/c_1  # Density in left half
    rho_2A  = gamma_2/c_2 # Density in right half
    bulk_1=bulk_1A
    bulk_2=bulk_2A
    rho_1=rho_1A
    rho_2=rho_2A
elif MaterialParams == 6:
    #print('Using Lurie's material parameters\n')
    gamma=1.0;
    gamma_1=gamma;
    gamma_2=gamma;
    c_1 = 0.55    # Sound speed (left)
    c_2 = 1.1  # Sound speed (right)
    bulk_1A = gamma_1*c_1  # Bulk modulus in left half
    bulk_2A= gamma_2*c_2  # Bulk modulus in right half
    rho_1A  = gamma_1/c_1  # Density in left half
    rho_2A  = gamma_2/c_2 # Density in right half
    bulk_1=bulk_1A
    bulk_2=bulk_2A
    rho_1=rho_1A
    rho_2=rho_2A
else:
    #This is know to work with the code that exists, so it is the default
    print('make sure to pick material params! Default used.\n')
    gamma=1.0;
    gamma_1=gamma;
    gamma_2=gamma+.1*gamma;
    c_1 = .6    # Sound speed (left)
    c_2 = 1.1  # Sound speed (right)
    bulk_1A = gamma_1*c_1  # Bulk modulus in left half
    bulk_2A= gamma_2*c_2  # Bulk modulus in right half
    rho_1A  = gamma_1/c_1  # Density in left half
    rho_2A  = gamma_2/c_2 # Density in right half
    bulk_1=bulk_1A
    bulk_2=bulk_2A
    rho_1=rho_1A
    rho_2=rho_2A


eps=.5;tau=.5;m=.5;n=.5;

alpha=.00001;beta= .00001;

timeInterfaceNum = 4;
t_0=0.0;t_F=timeInterfaceNum*tau;


#Splits regions up into checkerboard grid
def f_u(x,y,t,u1,u2):
    def M1(x,y): #These partition the function in space
        return (u1*(np.mod(x,eps)<m*eps) + u2*(np.mod(x,eps)>=m*eps));
    def M2(x,y):
        return (u2*(np.mod(x,eps)<m*eps) + u1*(np.mod(x,eps)>=m*eps));
    return (M1(x,y)*(np.mod(t,tau)<=n*tau)+M2(x,y)*(np.mod(t,tau)>n*tau)) #These partition the function in time

#Makes a linear interpolation between consecutive points and returns the interpolation value at xi
#def l(xi,xi1,xi2,y1,y2):
#    m=(y2-y1)/(xi2-xi1);
#    return y1+m*(xi-xi1)
    
#makes a 2d flat planar interpolation and returns the value at some point on the plane
#def p(xi,eta,xi1,xi2,eta1,eta2,y1,y2):
#    mz=(y2 - y1 )/(xi2 -xi1 );
#    mt=(y2 - y1 )/(eta2-eta1);
#    return y1+mz*(xi-xi1)+mt*(eta-eta1); 

#Sets up material properties to have linear gradients. We do not use that in this code though
#def py(x,z1,z2,t,t1,t2,u1,u2):
#    return (p(x,t,z1,z2,t1,t2,u1,u2)*((z1 <= x)*(x < z2)*(x-z1 <((z2-z1)/(t1-t2))*(t-t2)))+
#	    p(x,t,z2,z1,t2,t1,u1,u2)*((z1 <= x)*(x < z2)*(x-z1>=((z2-z1)/(t1-t2))*(t-t2))))
    


Wx=alpha;Wt=beta;
m1=m;n1=n;    



#Finds the total energy of an input state
#state is actually current_data where this function actually gets called, but the two objects have the same properties
def total_energy(state):
    totEnergy = 0.0;
    i = 0
    while i < (n_x):
        nxt = np.mod(i+1, n_x) #gets the correct next point in our space, which uses a wraparound boundary condition
#Energy is [((P_x)^2)/(1/rho) + ((V_x)^2)*K]/2 and K = rho*C so we may need to change the second term
#This energy term seems correct when considering the mathematical model of energy
        totEnergy += 0.5*((((state.q[0,nxt,0] - state.q[0,i,0])/SpaceStepSize)**2)/state.aux[0,i,0] + (((state.q[1,nxt,0] - state.q[1,i,0])/SpaceStepSize)**2)*(state.aux[1,i,0]**2)*state.aux[0,i,0])
#This is different from the derived term but appears to give a better result for energy
#        totEnergy += 0.5*((((state.q[0,nxt,0] - state.q[0,i,0])/SpaceStepSize)**2)/state.aux[0,i,0] + (((state.q[1,nxt,0] - state.q[1,i,0])/SpaceStepSize)**2)*state.aux[1,i,0])
        i += 1;
#need to update this function with numpy arrays when possible
#Try finding a pointer to these arrays instead of using indirection [0,i,0], get the vector before the while loop

    return totEnergy

#can implement this afterwords for fun
#def fourier():#Should use numpy fourier on the wave form to graph the spectrum


#We define all reflection and transmission coeficcients here
#reflection and transmission at spatial interface depends on which material that net flux is entering, and which material it is leaving
spatial_reflect_ab = ((gamma_1 - gamma_2)/(gamma_1 + gamma_2))**2
#spatial_transmit_ab = (gamma_2/gamma_1)*((2*gamma_1)/(gamma_1+gamma_2))**2
#We do not need the transmission coeficcients when calculating  energy exchanged at interfaces: energy that is not reflected is not exchanged
spatial_reflect_ba = ((gamma_2 - gamma_1)/(gamma_2 + gamma_1))**2
#spatial_transmit_ba = (gamma_1/gamma_2)*((2*gamma_2)/(gamma_2+gamma_1))**2

#Net temporal reflection and transmission depends on direction of net flux as we integrate over X
temporal_reflect = (0.5*((gamma_2/gamma_1) + gamma_1/gamma_2) - 1.0) #Same expression for material 2 to 1 or material 1 to 2
temporal_transmit = (0.5*((gamma_2/gamma_1) + gamma_1/gamma_2) + 1.0)

#Net amplification depends on current material and next material in each region
temporal_multiple_ab = 0.5*c_2/c_1
temporal_multiple_ba = 0.5*c_1/c_2

TimeToSpaceRatio = 0.8 #May need to be adjusted for various mesh sizes to get the right result, because Clawpack does not give direct access

#See notes about if statements
def energy_reflect(state):
    l = 0.0
    r = 0.0
    fluxlr = [l,r]
    flux = 0.0
    i = 0 
    nxt = np.mod(i+1, n_x)
    #Temporal interface from spatial laminate 1 to 2
    if ((np.mod(state.t,tau)<=n*tau) and (np.mod((state.t + SpaceStepSize*TimeToSpaceRatio),tau)>n*tau)) or ((np.mod(state.t,tau)>n*tau) and (np.mod((state.t + SpaceStepSize*TimeToSpaceRatio),tau)<=n*tau)):
        while i < (n_x): 
            if (state.aux[0,i,0] != state.aux[0,nxt,0]):
                #do averaging over spatial interface, then use that flux in the temporal switching calculation, with leftgoing and rightgoing fluxes partitioned into separate materials accordingly when doing the calculation
                flux = (((state.q[0,i,0] - state.q[0,i-1,0])*(state.q[1,i,0] - state.q[1,i-1,0]))/(SpaceStepSize**2))*(state.aux[1,i,0]**2)
                flux += (((state.q[0,nxt+1,0] - state.q[0,nxt,0])*(state.q[1,nxt+1,0] - state.q[1,nxt,0]))/(SpaceStepSize**2))*(state.aux[1,nxt,0]**2)
                flux = flux*0.5 



                if (state.aux[0,i,0] == rho_1): #Material 1 to 2 at interface
                    fluxlr[0] += flux*spatial_reflect_ab
                    fluxlr[1] -= flux*spatial_reflect_ab
                else:
                    fluxlr[0] += flux*spatial_reflect_ba
                    fluxlr[1] -= flux*spatial_reflect_ba



            else:#No need to average across the interface when calculating flux in region with constant material properties
                flux = (((state.q[0,nxt,0] - state.q[0,i,0])*(state.q[1,nxt,0] - state.q[1,i,0]))/(SpaceStepSize**2))*(state.aux[1,i,0]**2)
                if (state.aux[0,i,0] == rho_1): #Material 1 to 2 at interface
                        fluxlr[0] += flux*temporal_reflect*temporal_multiple_ab - (temporal_multiple_ba - 1.0)*flux*(flux < 0)
                        fluxlr[1] += (-1.0)*flux*temporal_reflect*temporal_multiple_ab + (temporal_multiple_ba - 1.0)*flux*(flux >= 0)
                        #transmitted energy multiplies: add new energy generated, avoiding double-counting of energy that was already there
                        #And  account for exchange by reflection separately    
                else: #Material 2 to 1 at interface
                        fluxlr[0] += flux*temporal_reflect*temporal_multiple_ba - (temporal_multiple_ba - 1.0)*flux*(flux < 0)
                        fluxlr[1] += (-1.0)*flux*temporal_reflect*temporal_multiple_ba + (temporal_multiple_ba - 1.0)*flux*(flux >= 0)
        i += 1
    #when not at a temporal interface, we just calculate energy exchanged at boundaries
    else:
        #we want the average of the fluxes on the two sides of the interfaces, not just the flux alone              
        while i < (n_x):
            print 'got to position 1' 
            print 'i is {0}'.format(i)
            if state.aux[0,i,0] != state.aux[0,nxt,0]:
                print 'got to position 2'
                #calculate flux
                flux = (((state.q[0,i,0] - state.q[0,i-1,0])*(state.q[1,i,0] - state.q[1,i-1,0]))/(SpaceStepSize**2))*(state.aux[1,i,0]**2)
                flux += (((state.q[0,nxt+1,0] - state.q[0,nxt,0])*(state.q[1,nxt+1,0] - state.q[1,nxt,0]))/(SpaceStepSize**2))*(state.aux[1,nxt,0]**2)
                flux = flux*0.5
                #partition flux 
                if (state.aux[0,i,0] == rho_1): #Material 1 to 2 at interface
                    fluxlr[0] += flux*spatial_reflect_ab #note that a negative flux will still be added and subtracted correctly
                    fluxlr[1] -= flux*spatial_reflect_ab
                    print 'got to position 3'
                else:
                    fluxlr[0] += flux*spatial_reflect_ba
                    fluxlr[1] -= flux*spatial_reflect_ba
                    print 'got to position 4'
        print 'i is {0}'.format(i)
        i += 1
#This is not quite correct: reflected energy must be subtracted from the group it was originally in, while transmitted energy is added, and then reflected energy should be added to the other group

#    def M1(x,y): #These partition the function in space
#        return (u1*(np.mod(x,eps)<m*eps) + u2*(np.mod(x,eps)>=m*eps));
#    def M2(x,y):
#        return (u2*(np.mod(x,eps)<m*eps) + u1*(np.mod(x,eps)>=m*eps));
#
    return fluxlr

#Calculates energy reflection and transmission using poynting vector
#Splits regions up into checkerboard grid
#def f_u(x,y,t,u1,u2):
#    def M1(x,y): #These partition the function in space
#        return (u1*(np.mod(x,eps)<m*eps) + u2*(np.mod(x,eps)>=m*eps));
#    def M2(x,y):
#        return (u2*(np.mod(x,eps)<m*eps) + u1*(np.mod(x,eps)>=m*eps));
#    return (M1(x,y)*(np.mod(t,tau)<=n*tau)+M2(x,y)*(np.mod(t,tau)>n*tau)) #These partition the function in time
        #This isn't quite the right way to do this yet
        #The poynting vector may take the form (1/(mu)(epsilon))*Phi_x*Psi_x
        #if BoundaryFluxFunction > 0
        #r_Energy += Transmissivity*BoundaryFluxFunction
        #r_Energy += Reflectivity*BoundaryFluxFunction
        #elif BoundaryFluxFunction < 0
        #r_Energy += Reflectivity*BoundaryFluxFunction
        #r_Energy += Transmissivity*BoundaryFluxFunction



#For use with the exponential limit curve, requires access to the global variables
base = np.amax([c_1,c_2])/np.amin([c_1,c_2])

def LimitCurve(Energy,T,T_0): #using initial energy, velocities and material geometry this will calculate the energy limit curve
    eLim = Energy*np.power(base,((T-T_0)/(n*tau)))
#np.exp(InitialEnergy, )
    return eLim




#Produces the initial conditions of the wave
def f_bump(z,z1,z2):
    def f4(z,z1,z2):
        return (z-z1)**2*(z-z2)**2
    return (f4(z,z1,z2)/f4((z1+z2)/2,z1,z2)*(z1<z)*(z<z2))


def setup(aux_time_dep=True,kernel_language='Fortran', use_petsc=False, outdir='./_output', 
          solver_type='classic', time_integrator='SSP104', lim_type=2, 
          disable_output=False, num_cells=(n_x, 1)):
    """
    Example python script for solving the 2d acoustics equations.
    """
    from clawpack import riemann

#    global Prevstep
    global l_Energy
    global r_Energy
    global InitEnergy

    if use_petsc:
        import clawpack.petclaw as pyclaw
    else:
        from clawpack import pyclaw

    if solver_type=='classic':
        solver=pyclaw.ClawSolver2D(riemann.vc_acoustics_2D)
        solver.dimensional_split=False
        solver.limiters = pyclaw.limiters.tvd.MC
    elif solver_type=='sharpclaw':
        solver=pyclaw.SharpClawSolver2D(riemann.vc_acoustics_2D)
        solver.time_integrator=time_integrator
        if time_integrator=='SSPLMMk2':
            solver.lmm_steps = 3
            solver.cfl_max = 0.25
            solver.cfl_desired = 0.24

    solver.bc_lower[0]=pyclaw.BC.periodic
    solver.bc_upper[0]=pyclaw.BC.periodic
    solver.bc_lower[1]=pyclaw.BC.periodic
    solver.bc_upper[1]=pyclaw.BC.periodic
    solver.aux_bc_lower[0]=pyclaw.BC.periodic
    solver.aux_bc_upper[0]=pyclaw.BC.periodic
    solver.aux_bc_lower[1]=pyclaw.BC.periodic
    solver.aux_bc_upper[1]=pyclaw.BC.periodic

    #uses parameters initialized at the start of the code to define space size and resolution
    x = pyclaw.Dimension(ax,bx,num_cells[0],name='x')
    y = pyclaw.Dimension(ay,by,num_cells[1],name='y')
    domain = pyclaw.Domain([x,y])

    num_eqn = 3
    num_aux = 2 # For the electromagnetic case the auxiliary variables are mu and epsilon
    state = pyclaw.State(domain,num_eqn,num_aux)

    grid = state.grid
    X, Y = grid.p_centers

#    
#    N_pl=1000
#    X_pl,T_pl = np.mgrid[ax:bx:1000*1j, t_0:t_F:1000*1j];
#    p_MG=plt.pcolor(X_pl,T_pl,f_u(X_pl,0.0,T_pl,c_1,c_2))
#    plt.savefig("_plots/Checkerboard.png",dpi=1000)


#    state.aux[0,:,:] = gamma/f_u(X,Y,0.0,c_1,c_2) # Density
##    state.aux[1,:,:] = f_u(X,Y,0.0,c_1  ,c_2  )    # Sound speed
#    state.aux[1,:,:] = f_u(X,Y,0.0,c_1  ,c_2)    # Sound speed

#sets initial material properties
    state.aux[0,:,:] = f_u(X,Y,0.0,rho_1,rho_2) # Density
#    state.aux[1,:,:] = f_u(X,Y,0.0,c_1  ,c_2  )    # Sound speed
    state.aux[1,:,:] = f_u(X,Y,0.0,c_1  ,c_2)    # Sound speed


    # Set initial condition
    x0 = -0.5; y0 = 0.
    r = np.sqrt((X-x0)**2 + (Y-y0)**2) # calculates a radial distance to a specified point (x0,y0)
    width = 0.10; rad = 0.25
    state.q[0,:,:] = f_bump(X,0.0,0.25) # sets the initial condition along the x direction
    state.q[1,:,:] = 0.
    state.q[2,:,:] = 0.
#    Prevstep = state.q
    #set up left energy and right energy here inside current_data
    InitEnergy = total_energy(state)
    l_Energy = InitEnergy*0.5
    r_Energy = InitEnergy*0.5

#!!Sets Local Material Properties State, outputs current wave state to buffer, calculates current energy and outputs it to CSV with current time step for plotting
    def DoBefore(solver,state):
        global l_Energy # Before each time step need to calculate reflection and transmission coefficients
        global r_Energy # 
        #TotEnergy = total_energy(state)
        flux = energy_reflect(state)
        l_Energy += flux[0] #This is not quite correct: reflected energy must be subtracted from the group it was originally in, while transmitted energy is added, and then reflected energy should be added to the other group
        r_Energy += flux[1]
        
#        state.aux[0,:,:] = f_u(X,Y,state.t,rho_1,rho_2)# Density
##        state.aux[1,:,:] = f_u(X,Y,state.t,c_1  ,c_2  ) # Sound speed       
#        state.aux[0,:,:] = gamma/f_u(X,Y,state.t,c_1  ,c_2  ) # Matching Impedances
        state.aux[0,:,:] = f_u(X,Y,state.t,rho_1,rho_2) # Density
        state.aux[1,:,:] = f_u(X,Y,state.t,c_1,c_2)    # Sound speed
        #print '{0},{1}'.format(state.t, TotEnergy) #original output for CSV, moved to setplot so that this doesn't have to be done ten thousand times
        #Prevstep = state.q   

    solver.before_step=DoBefore




    claw = pyclaw.Controller()
    claw.keep_copy = True
    if disable_output:
        claw.output_format = None
    claw.solution = pyclaw.Solution(state,domain)
    claw.solver = solver
    claw.outdir = outdir
    claw.tfinal = t_F
    claw.num_output_times = 100 #sets how many graphs are produced
    claw.write_aux_init = True
    claw.write_aux=True
    claw.setplot = setplot
    if use_petsc:
        claw.output_options = {'format':'binary'}

    return claw


def setplot(plotdata):#maybe it is possible to add state here and use it instead of current_data? Not sure if that would break the calling sequence
    """ 
    Plot solution using VisClaw.

    This example shows how to mark an internal boundary on a 2D plot.
    """ 

    global l_Energy
    global r_Energy
    global InitEnergy

    from clawpack.visclaw import colormaps
#Plot all new functions in setplot
    plotdata.clearfigures()  # clear any old figures,axes,items data

#The commented section here are figures from the example code that we probably don't need    
    # Figure for pressure
#   plotfigure = plotdata.new_plotfigure(name='Pressure', figno=0)

    # Set up for axes in this figure:
#    plotaxes = plotfigure.new_plotaxes()
#    plotaxes.title = 'Pressure'
#    plotaxes.scaled = True      # so aspect ratio is 1
#    plotaxes.afteraxes = mark_interface
    #plotaxes.afteraxes = other_function #Can insert functions to plot this way
#    plotaxes.ylimits=[ay,by]
#    plotaxes.xlimits=[ax,by] #This is a hack because bx doesn't work here for some reason although they should be the same value
    #plotaxes.xlimits=[ax,bx]

    # Set up for item on these axes:
#    plotitem = plotaxes.new_plotitem(plot_type='2d_pcolor')
#    plotitem.plot_var = 0
#    plotitem.pcolor_cmap = colormaps.yellow_red_blue
#    plotitem.add_colorbar = True
#    plotitem.pcolor_cmin = 0.0
#    plotitem.pcolor_cmax=1.0
    
    # Figure for x-velocity plot
#    plotfigure = plotdata.new_plotfigure(name='x-Velocity', figno=1)

    # Set up for axes in this figure:
#    plotaxes = plotfigure.new_plotaxes()
#    plotaxes.title = 'u'
#    plotaxes.afteraxes = mark_interface

#    plotitem = plotaxes.new_plotitem(plot_type='2d_pcolor')
#    plotitem.plot_var = 1
#    plotitem.pcolor_cmap = colormaps.yellow_red_blue
#    plotitem.add_colorbar = True
#    plotitem.pcolor_cmin = -0.3
#    plotitem.pcolor_cmax=   0.3


    # Would be nice to graph phi, psi, E, H and Energy Density as functions of the propagating wave
    plotfigure = plotdata.new_plotfigure(name='Psi', figno=0)

    # Set up for axes in this figure:
    plotaxes = plotfigure.new_plotaxes()
    plotaxes.title = 'Potential'
    plotaxes.scaled = True      # so aspect ratio is 1
    plotaxes.afteraxes = mark_interface
    plotaxes.ylimits = [ay,by]
    plotaxes.xlimits = [ax,bx]

    # Set up for item on these axes:
    plotitem = plotaxes.new_plotitem(plot_type='1d_from_2d_data')
    plotitem.plot_var = 0

    def q_y0(current_data):
        x=current_data.x
        qy=current_data.var
        return x,qy
    plotitem.map_2d_to_1d = q_y0

    def add_plot(current_data):
        global l_Energy
        global r_Energy
        global InitEnergy

        x = current_data.x;
        y = current_data.y;
        a=current_data.aux[1]
        
        print current_data.t
        Ones=np.ones(x.shape)
        plt.plot(x,.5*Ones, 'k')
        #T1plt.plot(x,gamma/f_u(x,y,current_data.t,rho_1,rho_2),'-r')
        plt.plot(x,      f_u(x,y,current_data.t,c_1,c_2),'-r')
        plt.plot(x,a,'-g')
        plt.plot(x, ) #What is going on with this?
        plt.title('EM Potentials and Waves')
        #plt.plot(x,gamma/f_u(x,y,current_data.t,c_1,c_2),'-g')
        if current_data.t == 1.475: #Magic numbers for m=n=0.5 and 10000 points in x:need to think about how to generalize
            time = 0.0
            Energy = total_energy(current_data)
            while time <= 2.5:
                print 'LimitCurve {0}'.format(LimitCurve(Energy,time,current_data.t))
                time += 0.025
        #These print statements are meant for postprocessing by a perl file to produce a CSV which represents energy over time
        print 'EnergyOutput {0},{1},{2},{3},{4}'.format(current_data.t, total_energy(current_data), l_Energy, r_Energy, (l_Energy+r_Energy))
    plotaxes.afteraxes = add_plot #Allows us to plot additional functions of current_data
    #Next we output values for energy for use in a CSV to be graphed in another program
    plotaxes.xlimits=[ax,bx]   
    plotaxes.ylimits=[0.0,1.8] #should really be the max of wave amplitude, though that may change depending on the energy

#Not a full plotting example: use the others above which are also commented out to figure out the right way to do this
#Would be nice to also plot the fourier transform of state.q to allow visual inspection of information loss as the wave approaches the resolution limit of the mesh
    # Fourier transform over time of the wave
#    plotfigure = plotdata.new_plotfigure(name='Frequency Spectrum', figno=4)

    # Set up for axes in this figure:
#    plotaxes = plotfigure.new_plotaxes()
#    plotaxes.title = 'Frequency Spectrum'
#    plotaxes.scaled = True      # so aspect ratio is 1
#    plotaxes.ylimits=[ay,by] #What axes do I want in my frequency domain?
    #Half the nyquist limit probably
#    plotaxes.xlimits=[ax,by] #This shouldn't work but it does
    


    # Parameters used only when creating html and/or latex hardcopy
    # e.g., via visclaw.frametools.printframes:

    plotdata.printfigs = True                # print figures
    plotdata.print_format = 'png'            # file format
    plotdata.print_framenos = 'all'          # list of frames to print
    plotdata.print_fignos = 'all'            # list of figures to print
    plotdata.html = True                     # create html files of plots?
    plotdata.html_homelink = '../README.html'   # pointer for top of index
    plotdata.html_movie = 'JSAnimation'      # new style, or "4.x" for old style
    plotdata.latex = True                    # create latex file of plots?
    plotdata.latex_figsperline = 2           # layout of plots
    plotdata.latex_framesperline = 1         # layout of plots
    plotdata.latex_makepdf = False           # also run pdflatex?
    
    return plotdata

def mark_interface(current_data):
    import matplotlib.pyplot as plt
    plt.plot((0.,0.),(-1.,1.),'-k',linewidth=2)


if __name__=="__main__":
    from clawpack.pyclaw.util import run_app_from_main
    output = run_app_from_main(setup,setplot)
    #so the way setplot gets called when the program runs is determined in run_app_from_main
