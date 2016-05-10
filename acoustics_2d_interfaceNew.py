#!/usr/bin/env python
# encoding: utf-8
r"""
1D acoustics equation with variable coeficcients in space and time
==============================================

Solve the variable-coefficient acoustics equations in 2D:

.. math:: 
    p_t + K(x,t) (u_x) & = 0 \\ 
    u_t + p_x / \rho(x,t) & = 0 \\

Here p is the pressure, (u,v) is the velocity, :math:`K(x,y)` is the bulk modulus,
and :math:`\rho(x,y)` is the density.

This example shows how to solve a problem with variable coefficients.

Note that it is equivalent to maxwell's system in 1D, if we 
"""
 
import numpy as np
import matplotlib.pyplot as plt

#here we define the size of the simulated space
ax=0.0;bx=1.0;n_x=500
SpaceStepSize = (bx - ax)/n_x #note that this is only accurate for an evenly-spaced grid
ay=0.0;by=1.0;
#n_y is just 1 because we have a 1D system, not a 2D system


MaterialParams = 1; #set this variable to pick which material properties we have

if MaterialParams == 1:
    print('default with slight reflection\n') #This is known to work with the code that exists
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
    print('No impedance mismatch or reflection\n')
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
    print('Large reflection\n')
    gamma=1.0;
    gamma_1=gamma;
    gamma_2=gamma+.2*gamma;
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
    print('Large wave speed ratio and reflection\n')
    gamma=1.0;
    gamma_1=gamma;
    gamma_2=gamma+.2*gamma;
    c_1 = .2    # Sound speed (left)
    c_2 = 1.1  # Sound speed (right)
    bulk_1A = gamma_1*c_1  # Bulk modulus in left half
    bulk_2A= gamma_2*c_2  # Bulk modulus in right half
    rho_1A  = gamma_1/c_1  # Density in left half
    rho_2A  = gamma_2/c_2 # Density in right half
    bulk_1=bulk_1A
    bulk_2=bulk_2A
    rho_1=rho_1A
    rho_2=rho_2A
elif MaterialParams == 5:
    print('sanity check: no reflection no wavespeed so no energy change\n')
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
else: #This is know to work with the code that exists, so it is the default
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




def f_u(x,y,t,u1,u2):
    def M1(x,y):
        return (u1*(np.mod(x,eps)<m*eps) + u2*(np.mod(x,eps)>=m*eps));
    def M2(x,y):
        return (u2*(np.mod(x,eps)<m*eps) + u1*(np.mod(x,eps)>=m*eps));
    return (M1(x,y)*(np.mod(t,tau)<=n*tau)+M2(x,y)*(np.mod(t,tau)>n*tau))

def l(xi,xi1,xi2,y1,y2):
    m=(y2-y1)/(xi2-xi1);
    return y1+m*(xi-xi1)
    
def p(xi,eta,xi1,xi2,eta1,eta2,y1,y2):
    mz=(y2 - y1 )/(xi2 -xi1 );
    mt=(y2 - y1 )/(eta2-eta1);
    return y1+mz*(xi-xi1)+mt*(eta-eta1); 
    
def py(x,z1,z2,t,t1,t2,u1,u2):
    return (p(x,t,z1,z2,t1,t2,u1,u2)*((z1 <= x)*(x < z2)*(x-z1 <((z2-z1)/(t1-t2))*(t-t2)))+
	    p(x,t,z2,z1,t2,t1,u1,u2)*((z1 <= x)*(x < z2)*(x-z1>=((z2-z1)/(t1-t2))*(t-t2))))
    



Wx=alpha;Wt=beta;

m1=m;n1=n;    





def f_bump(z,z1,z2):
    def f4(z,z1,z2):
        return (z-z1)**2*(z-z2)**2
    return f4(z,z1,z2)/f4((z1+z2)/2,z1,z2)*(z1<z)*(z<z2)



def setup(aux_time_dep=True,kernel_language='Fortran', use_petsc=False, outdir='./_output', 
          solver_type='classic', time_integrator='SSP104', lim_type=2, 
          disable_output=False, num_cells=(n_x, 1)):
    """
    Example python script for solving the 2d acoustics equations.
    """
    from clawpack import riemann

    global Prevstep

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

    x = pyclaw.Dimension(ax,bx,num_cells[0],name='x')
    y = pyclaw.Dimension(ay,by,num_cells[1],name='y')
    domain = pyclaw.Domain([x,y])

    num_eqn = 3
    num_aux = 2 # density, sound speed
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

    state.aux[0,:,:] = f_u(X,Y,0.0,rho_1,rho_2) # Density
#    state.aux[1,:,:] = f_u(X,Y,0.0,c_1  ,c_2  )    # Sound speed
    state.aux[1,:,:] = f_u(X,Y,0.0,c_1  ,c_2)    # Sound speed


    # Set initial condition
    x0 = -0.5; y0 = 0.
    r = np.sqrt((X-x0)**2 + (Y-y0)**2)
    width = 0.1; rad = 0.25
    state.q[0,:,:] = f_bump(X,.25,.75)
    state.q[1,:,:] = 0.
    state.q[2,:,:] = 0.
    Prevstep = state.q

#!!Sets Local Material Properties State, outputs current wave state to buffer, calculates current energy and outputs it to CSV with current time step for plotting
    def DoBefore(solver,state):
        #global Prevstep
        TotEnergy = 0.0;
#        state.aux[0,:,:] = f_u(X,Y,state.t,rho_1,rho_2);# Density
##        state.aux[1,:,:] = f_u(X,Y,state.t,c_1  ,c_2  ); # Sound speed       
#        state.aux[0,:,:] = gamma/f_u(X,Y,state.t,c_1  ,c_2  ); # Matching Impedances
        i = 0;
        while i < (n_x-1):
           TotEnergy += (state.aux[0,i,0]*((state.q[0,i+1,0] - state.q[0,i,0])/SpaceStepSize)**2) + ((((state.q[1,i+1,0] - state.q[1,i,0])/SpaceStepSize)**2)/state.aux[1,i,0])
           i += 1;
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
    claw.num_output_times = 100
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

    plotdata.clearfigures()  # clear any old figures,axes,items data
    
    # Figure for pressure
    plotfigure = plotdata.new_plotfigure(name='Pressure', figno=0)

    # Set up for axes in this figure:
    plotaxes = plotfigure.new_plotaxes()
    plotaxes.title = 'Pressure'
    plotaxes.scaled = True      # so aspect ratio is 1
    plotaxes.afteraxes = mark_interface
    plotaxes.ylimits=[ay,by]
    plotaxes.xlimits=[ax,by]

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
   



    plotfigure = plotdata.new_plotfigure(name='1D-Pressure', figno=3)

    # Set up for axes in this figure:
    plotaxes = plotfigure.new_plotaxes()
    plotaxes.title = 'Pressure'
    plotaxes.scaled = True      # so aspect ratio is 1
    plotaxes.afteraxes = mark_interface
    plotaxes.ylimits=[ay,by]
    plotaxes.xlimits=[ax,bx]

    # Set up for item on these axes:
    plotitem = plotaxes.new_plotitem(plot_type='1d_from_2d_data')
    plotitem.plot_var = 0

    def q_y0(current_data):
        x=current_data.x
        qy=current_data.var
        return x,qy
    plotitem.map_2d_to_1d=q_y0



    def add_plot(current_data):
        x = current_data.x;y=current_data.y;a=current_data.aux[1]
        
        print current_data.t
 	Ones=np.ones(x.shape)
        plt.plot(x,.5*Ones, 'k')
        #T1plt.plot(x,gamma/f_u(x,y,current_data.t,rho_1,rho_2),'-r')
        plt.plot(x,      f_u(x,y,current_data.t,c_1,c_2),'-r')
        plt.plot(x,a,'-g')
        plt.title('Stuff')
        #plt.plot(x,gamma/f_u(x,y,current_data.t,c_1,c_2),'-g')
    plotaxes.afteraxes = add_plot
    
    plotaxes.xlimits=[ax,bx]   
    plotaxes.ylimits=[0.0,1.3] 



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
