# MQP
Wave Equation Solvers for Senior Thesis in Python

Look at the makefile to see examples of what command line prompts are needed to
 run the code.

The output to the command line is a CSV of current time and current total energy
 of the wave, but Clawpack itself has print statements internally that I don't
want to mess with, so we have Postprocess.pl to remove those statements

##Dependencies:
The clawpack documentation for getting things running is excellent
[read it here](http://www.clawpack.org/installing.html#installation-instructions)
Especially read the part about setting the environment variables: put those
command line prompts in .bashrc in your home folder.

But my code may have additional dependencies, so I have listed everything I
remember here:
1. Clawpack-5.3.1
2. gfortran or an equivalent Fortran compiler
3. Python 2.7
4. Perl 5.22.2
5. Various Python math libraries discussed in the Clawpack documentation

Also if using emclaw for more general EM Simulation you will need
1. petsc
2. petsc4py
3. clawpack-4.6.3 instead of 5.3.1
