# MQP
Wave Equation Solvers for Senior Thesis in Python

## Current Workflow
1. During "make graphs" the output to the command line is saved as a text file post-processed by perl and keeps track of of current time, current total energy of the wave, and the energy in leftgoing and rightgoing families of characteristics.
2. A separate CSV from perl post-processing contains a limit curve for the exponential growth.
3. For now to change the material properties you need to set a switch statement in the code.
4. Once you get EnergyVals.csv and LimitCurve.csv, import to excel or an equivalent spreadsheet program to graph and visualize it.
5. The actual wave propagation video is stored in the output folder Plots. Other things can be graphed this way as well. It can be accessed by the URL printed to the command line on completion of 'make graphs'

##Dependencies:

To install most of the dependencies, we can use
pip install -r requirements.txt
in the terminal or check out the dependencies in requirements.txt if you want to install another way, in case pip fails. As of 2018 it seems pip works in virtualenv for all of the requirements except gfortran and perl. The dependencies should be installed in virtualenv for portability, and some instructions for that are below. pip 8.x and Python 2.7.12 must already be installed if you want to run things locally outside virtualenv, otherwise the version numbers are not as important. All these programs can most likely be installed by apt-get too, or at least have their own documentation.
To set up virtualenv in a particular location with a particular python, use
virtualenv --python=/usr/bin/python2.7 (/path to this project/ven)
This lets you select from a number of installed python versions. The name of the virtual environmet folder is assumed to be ven by gitignore so might as well use it.


The clawpack documentation for getting things running is excellent
[read it here](http://www.clawpack.org/installing.html#installation-instructions)
Especially read the part about setting the environment variables: put those
command line prompts in .bashrc in your home folder. The environment variables CLAW and PYTHONPATH allow clawpack and python to be imported in those respective folders
It should look something like this:
export CLAW=/home/michael/GitFolder/MQP/clawpack-5.4.0
export FC=gfortran
export PYTHONPATH=/home/michael/GitFolder/MQP/clawpack-5.4.0/clawpack

pythonpath should look like this because it must be imported inside clawpack

More info about pythonpath:
https://stackoverflow.com/questions/15109548/set-pythonpath-before-import-statements


The default system python path can be found by typing 'which python', but if you need 2.7.12 and don't have it, you need to use virtualenv, so you install it separately and point pythonpath to the folder this project is in, and then set up pythonpath to point to the directory of the new virtual environment. you can do this with
virtualenv --python=/usr/bin/python2.7 (this arg is your pwd in this folder, where you should point pythonpath)

For portability, we use virtualenv (with python 2.7.12) to run the clawpack code, which makes a virtual python environment in the current directory, so that this project doesn't interfere with the default Python version of the computer. To set up virtualenv, see the instructions further below.




Other things to note:
1. gfortran can be replaced with an equivalent fortran compiler. As of 2018 it seems gfortran is easier to get with apt-get instead of pip
2. Make sure pip is installed and up to date so that everything installs properly - more precisely installation works best with pip
3. in Virtualenv, Python 2.7.12 was used for the project. Other Python versions may be compatible, but it is recommended to use the one the project originally had.
4. You might need to install Clawpack 5.4  and gfortran via apt get instead of pip
5. You might need to do sudo pip install clawpack after installing python for things to work
6. I should probably make all this into a script later on
7. On one computer I installed to, six.py was missing for some reason. It is forseeable that other people's computers may be missing dependencies they ought to have by default
8. on Ubuntu's default Python (in my case 2.7.12) tkinter might not be working, which might be solved via removing the virtual environment, sudo apt-get install tk-dev python-tk, and then adding the virtual environment back so it makes python with tkinter. If not using virtualenv, or if that doesn't work, you may need to compile python with tkinter manually unfortunately.


To set the Python version in virtualenv, see the following conversation online:
http://stackoverflow.com/questions/1534210/use-different-python-version-with-virtualenv
or see the '-p' arguement in the man page for virtualenv

## Using virtualenv
In general, to create a virtual environment use
1. virtualenv name

Note .gitignore is set to ignore a virtual environment named ven. Activate virtual environment with
2. . ven/bin/activate

deactivate virtual environment with
3. deactivate

install dependencies for the project with
4. pip install -r requirements.txt

When virtualenv is activated, the makefile can run the code with
5. make graphs

Important: for python, tabs need to be converted into four spaces, but for the makefile they need to be actual tabs. if the files are called 'makefile' vim will recognize them and change behavior. Just in case, see the provided example for your .vimrc file that you should create in your user directory. If you need tabs in other circumstances, ctrl-V-tab will produce an actual tab regardless of vimrc.

Obviously vim is not needed for running things but I like it.
