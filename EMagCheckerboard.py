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
ax=0.0;bx=1.0;n_x=200
#ax=0.0;bx=1.0;n_x=5000 #higher resolution allows us to avoid non-physical resolution limit of power amplification for longer. Mesh refinement would be better.
SpaceStepSize = (bx - ax)/n_x #note that this is only accurate for an evenly-spaced grid. For mesh refinement an alternative must be found
ay=0.0;by=1.0;n_y=1
#n_y is just 1 because we have a 1D system, not a 2D system, though the original code was for 2D


MaterialParams = 4 #set this variable to pick which material properties we have
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

t_0=0.0;t_F=5*tau;



#Splits regions up into checkerboard grid
def f_u(x,y,t,u1,u2):
    def M1(x,y): #These partition the function in space
        return (u1*(np.mod(x,eps)<m*eps) + u2*(np.mod(x,eps)>=m*eps));
    def M2(x,y):
        return (u2*(np.mod(x,eps)<m*eps) + u1*(np.mod(x,eps)>=m*eps));
    return (M1(x,y)*(np.mod(t,tau)<=n*tau)+M2(x,y)*(np.mod(t,tau)>n*tau)) #These partition the function in time

#Makes a linear interpolation between consecutive points and returns the interpolation value at xi
def l(xi,xi1,xi2,y1,y2):
    m=(y2-y1)/(xi2-xi1);
    return y1+m*(xi-xi1)
    
#makes a 2d flat planar interpolation and returns the value at some point on the plane
def p(xi,eta,xi1,xi2,eta1,eta2,y1,y2):
    mz=(y2 - y1 )/(xi2 -xi1 );
    mt=(y2 - y1 )/(eta2-eta1);
    return y1+mz*(xi-xi1)+mt*(eta-eta1); 
    
def py(x,z1,z2,t,t1,t2,u1,u2):
    return (p(x,t,z1,z2,t1,t2,u1,u2)*((z1 <= x)*(x < z2)*(x-z1 <((z2-z1)/(t1-t2))*(t-t2)))+
	    p(x,t,z2,z1,t2,t1,u1,u2)*((z1 <= x)*(x < z2)*(x-z1>=((z2-z1)/(t1-t2))*(t-t2))))
    


Wx=alpha;Wt=beta;
m1=m;n1=n;    



#Finds the total energy of an input state
def total_energy(state):
    totEnergy = 0.0;
    i = 0
    while i < (n_x):
        nxt = np.mod(i+1, n_x) #gets the correct next point in our space, which uses a wraparound boundary condition
        totEnergy += (((state.q[0,nxt,0] - state.q[0,i,0])/SpaceStepSize)**2)/state.aux[0,i,0] + (((state.q[1,nxt,0] - state.q[1,i,0])/SpaceStepSize)**2)*state.aux[1,i,0]
        i += 1;
        print '{0},{1}'.format(state.t, totEnergy) #got rid of the other print statement for debugging purposes. Need to modify the perl postprocess accordingly       
    return totEnergy

#def fourier():#Should use numpy fourier on the wave form to graph the spectrum

#def energy_reflect(): #Calculates energy reflection at the spatial boundaries using average, weighted average, and extrapolation of the poynting vector for estimates of the poynting vector at 

#def energy_t_reflect(): # Calculates energy reflection at the step before a temporal boundary, including the effect of the boundary in terms of amplification

#def LimitCurve(InitialEnergy,M,N,V1,V2,T): #using initial energy, velocities and material geometry this will calculate the limit curve for energy at time T
        #
#    return 0



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
        #global Prevstep
        global l_Energy #Before each time step need to calculate reflection and transmission coefficients
        global r_Energy
        TotEnergy = 0.0;
        #TotEnergy = total_energy(state);

#        state.aux[0,:,:] = f_u(X,Y,state.t,rho_1,rho_2);# Density
##        state.aux[1,:,:] = f_u(X,Y,state.t,c_1  ,c_2  ); # Sound speed       
#        state.aux[0,:,:] = gamma/f_u(X,Y,state.t,c_1  ,c_2  ); # Matching Impedances
        
        i = 0;
        while i < (n_x):
           nxt = np.mod(i+1, n_x) #We have a wraparound boundary condition so we must calculate energy around the whole loop
           #It is worth it to store the arrays involved in numpy arrays to avoid having to loop in python
           #Also this really ought to be a separate function that is called in do_before

           #Try finding a pointer to these arrays instead of using indirection [0,i,0], get the vector before the while loop
           TotEnergy += (((state.q[0,nxt,0] - state.q[0,i,0])/SpaceStepSize)**2)/state.aux[0,i,0] + (((state.q[1,nxt,0] - state.q[1,i,0])/SpaceStepSize)**2)*state.aux[1,i,0]
           i += 1; #make these more eficcient with numpy arrays
           #insert an array to save the previous time step/state.q, state.h, then put setplot in here and turn it off outside - go see why it only works 100 times instead of every time step in the other code also.
           #We already calculate the derivatives of the potentials above, so why not use those to graph D and B instead of just the potentials here
        state.aux[0,:,:] = f_u(X,Y,state.t,rho_1,rho_2) # Density
        state.aux[1,:,:] = f_u(X,Y,state.t,c_1,c_2)    # Sound speed
        #print '{0},{1},{2}'.format(state.t, TotEnergy, n_x) #for debugging
        print '{0},{1}'.format(state.t, TotEnergy)
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


def setplot(plotdata):
    """ 
    Plot solution using VisClaw.

    This example shows how to mark an internal boundary on a 2D plot.
    """ 

    from clawpack.visclaw import colormaps
#Plot all new functions in setplot
    plotdata.clearfigures()  # clear any old figures,axes,items data
    
    # Figure for pressure
    plotfigure = plotdata.new_plotfigure(name='Pressure', figno=0)

    # Set up for axes in this figure:
    plotaxes = plotfigure.new_plotaxes()
    plotaxes.title = 'Pressure'
    plotaxes.scaled = True      # so aspect ratio is 1
    plotaxes.afteraxes = mark_interface
    #plotaxes.afteraxes = other_function #Can insert functions to plot this way
    plotaxes.ylimits=[ay,by]
    plotaxes.xlimits=[ax,by] #This is a hack because bx doesn't work here for some reason although they should be the same value
    #plotaxes.xlimits=[ax,bx]

    # Set up for item on these axes:
    plotitem = plotaxes.new_plotitem(plot_type='2d_pcolor')
    plotitem.plot_var = 0
    plotitem.pcolor_cmap = colormaps.yellow_red_blue
    plotitem.add_colorbar = True
    plotitem.pcolor_cmin = 0.0
    plotitem.pcolor_cmax=1.0
    
    # Figure for x-velocity plot
    plotfigure = plotdata.new_plotfigure(name='x-Velocity', figno=1)

    # Set up for axes in this figure:
    plotaxes = plotfigure.new_plotaxes()
    plotaxes.title = 'u'
    plotaxes.afteraxes = mark_interface

    plotitem = plotaxes.new_plotitem(plot_type='2d_pcolor')
    plotitem.plot_var = 1
    plotitem.pcolor_cmap = colormaps.yellow_red_blue
    plotitem.add_colorbar = True
    plotitem.pcolor_cmin = -0.3
    plotitem.pcolor_cmax=   0.3


    plotfigure = plotdata.new_plotfigure(name='Psi', figno=3)

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
        #output energy, limit curve, l_Energy and r_Energy to a CSV here
    plotaxes.afteraxes = add_plot
    
    plotaxes.xlimits=[ax,bx]   
    plotaxes.ylimits=[0.0,1.8] #should really be the max of wave amplitude, though that may change depending on the 
 
    # Fourier transform over time of the wave
    plotfigure = plotdata.new_plotfigure(name='Frequency Spectrum', figno=4)

    # Set up for axes in this figure:
    plotaxes = plotfigure.new_plotaxes()
    plotaxes.title = 'Frequency Spectrum'
    plotaxes.scaled = True      # so aspect ratio is 1
    plotaxes.ylimits=[ay,by] #What axes do I want in my frequency domain?
    #Half the nyquist limit probably
    plotaxes.xlimits=[ax,by] #This shouldn't work but it does
    


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
