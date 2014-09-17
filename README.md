Schroe.py v0.5
==============

https://github.com/heedmane/schroepy/
Licensed under GPLv2 (http://www.gnu.org/licenses/gpl-2.0.html)

Python port of the Mathematica script specified in arXiv:hep-ph/9811453

Author:
Hector E. Martinez, 
Physik-Department T30f,
TU Muenchen,
Garhing, Germany.
hector.martinez@tum.de

cfunctions.c include code written by Thomas Rosenhammer
SClib licensed under GPLv2 https://github.com/drestebon/sclib/


Requirements
------------
*nix OS (Linux, Mac OS X)

For basic pure python usage (available in the PPY/ folder):
 * Python
 * IPython
 * NumPy
 * SciPy
 * Matplotlib

To use SClib (available in the SC/ folder):
 * GCC
 * Make

I haven't tested the script in windows but if you can install the programs in the list
in windows it should run with no problems.


Installation 
------------

* Download (or clone) the current/ folder to your machine.

* If you want to use the SC version inside the SC/src/ foler rename the Makefile correspondiing to your OS:

    cp Makefile_mac Makefile
or 
    cp Makefile_linux Makefile

* Inside PPy/src/ or SC/src/ folder start IPython. It is convenient to use the option 'inline' of matplotlib to see the plots. If you run the qtconsole run with
 
    ipython qtconsole --pylab="inline"

or inside a IPyhton notebook use the magic

    %matplotlib inline

* Within IPython type 'run SChroe.py' This will run the script with the example potential.


Basic Usage 
-----------
* Modify the SChroe.py and/or potential.h to define your potential, see the comments therein. 

* To calculate and store the eigenvalue and the wavefunction:

    Enl,ynl = solve_schrodingerPy(elow,eup,n,l)

or

    Enl,ynl = solve_schrodingerC(elow,eup,n,l)


* To only compute and store the nth (n=0,1,2,...) eigenvalue 

    Enl = eigenvalueC(elow,eup,n,l)
or
    Enl = eigenvaluePy(elow,eup,n,l)

The fist function is faster because uses SClib, the second one is the Python version of the same algorithm.

* To calculate and store the reduced wavefunction corresponding to this eigenvalue:
    ynl = wavefunctionC(Enl,l) 
or  
    ynl = wavefunctionPy(Enl,l)

* Check orthonormality:
    1-wfint(wfpro2(ynl,ynl)) 
and 
    wfint(wfpro2(ynl,yml)) 
where m != n.

* plotting: the function 
    wfplot2([ynl],[xl,xf]) 
plots the reduced wavefunction between xl and xf, or you can use

    wfplot2([ynl],-1) 
to plot it in its whole range. Also you can use

    wfplot2([ynl,yml,...],-1)

to plot two or more wavefunctions.

* save wavefunctions to disk with

    np.save('name',ynl) 

load with 
    ynl = np.load('name')

* See the notebook in the example0 folder.

* If you want the potential to depend on more parameters than r see the notebook inside the example2 folder and the comments therein.

* If you make changes in potential.h or cfunctions.c you may need to erase cfunctions.so

