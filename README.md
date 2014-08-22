---------------
---------------
Schroe.py v1.0
---------------
---------------

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
------------
Unix-like OS (Linux, Mac OS X)

For basic usage:
 - Python
 - IPython
 - NumPy
 - SciPy
 - Matplotlib

To use SClib:
 - GCC
 - Make

Recommended:
- IPython qtconsole

I haven't tested the script in windows but if you can install the programs in the list
in windows it should run with no problems.

** If you do not want to use SClib erase that part from the script file. **

Installation 
------------
------------

- Download (or clone) to your machine.

- Rename the Makefile, for instance in Mac OS X using gcc4.2 you need rename 'Makefile_mac' to Makefile, for Linux use 'Makefile_linux'

- Inside src/ start IPython with the option --pylab="inline" (We recommend to run the script in the IPython qtconsole)
 
    ' ipython qtconsole --pylab="inline" '

- Within IPython type 'run SChroe.py' This will run the script with the example potential.


Basic Usage 
-----------
-----------
- Modify the SChroe.py and potential.h to define your potential, see the comments therein. 
 
- To compute and store the nth (n=0,1,2,...) eigenvalue 'Enl = eigenvalueC(elow,eup,n,l)' or 'Enl = eigenvaluePy(elow,eup,n,l)'
- The fist function is faster because uses SClib, the second one is the Python version of the same algorithm.

- To calculate and store the reduced wavefunction corresponding to this eigenvalue ynl = wavefunction(Enl,l)

- Check orthonormality 1-wfint(wfpro2(ynl,ynl)) or wfint(wfpro2(ynl,yml)) where m != n.

- plotting: wfplot2([ynl],[xl,xf]) plots the reduced wavefunction between xl and xf, or wfplot2([ynl],-1) to plot it in its whole range,  wfplot2([ynl,yml],-1) to plot two wavefunctions.

- save wavefunctions to disk with np.save('name',ynl) load with ynl = np.load('name')

- If you want the potential to depend on more parameters than r see example1 and the comments therein.

- If you make changes in potential.h or cfunctions.c you may need to erase cfunctions.so


