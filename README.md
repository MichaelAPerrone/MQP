# MQP
Wave Equation Solvers for Senior Thesis in Python

## Current Workflow
1. During "make graphs" the output to the command line is saved as a text file post-processed by perl and keeps track of of current time, current total energy of the wave, and the energy in leftgoing and rightgoing families of characteristics.
2. A separate CSV from perl post-processing contains a limit curve for the exponential growth.
3. For now to change the material properties you need to set a switch statement in the code.
4. Once you get EnergyVals.csv and LimitCurve.csv, import to excel or an equivalent spreadsheet program to graph and visualize it.
5. The actual wave propagation video is stored in the output folder Plots. Other things can be graphed this way as well. It can be accessed by the URL printed to the command line on completion of 'make graphs'

##Dependencies:
The clawpack documentation for getting things running is excellent
[read it here](http://www.clawpack.org/installing.html#installation-instructions)
Especially read the part about setting the environment variables: put those
command line prompts in .bashrc in your home folder. There is a py

For portability, we use virtualenv to run the clawpack code, which makes a virtual python environment in the current directory, so that this project doesn't interfere with the default Python version of the computer.

Once Virtualenv is active, we can install most of the dependencies for the project with
pip install -r requirements.txt
or check out what the dependencies are by opening requirements.txt if you prefer to install another way.

Other things to note:
1. gfortran can be replaced with an equivalent fortran compiler
2. Make sure pip is installed and up to date so that everything installs properly - more precisely installation works best with pip
3. in Virtualenv, Python 2.7.12 was used for the project. Other Python versions may be compatible, but it is recommended to use the one the project originally had.
4. You might need to install Clawpack 5.4  and gfortran via apt get instead of pip
5. You might need to do sudo pip install clawpack after installing python for things to work
6. I should probably make all this into a script later on

To set the Python version in virtualenv, see the following conversation online:
http://stackoverflow.com/questions/1534210/use-different-python-version-with-virtualenv
or see the '-p' arguement in the man page for virtualenv

## Using virtualenv
In general, to create a virtual environment use
1. virtualenv name

activate virtual environment with
2. . ven/bin/activate

deactivate virtual environment with
3. deactivate

install dependencies for the project with
4. pip install -r requirements.txt

When virtualenv is activated, the makefile can run the code with
5. make graphs

Note to self: for python, tabs need to be converted into four spaces, but for the makefile they need to be actual tabs. In Vim, shift tab gives you tab, and 'set et' and 'set sw=4' makes tab into 4 spaces.
