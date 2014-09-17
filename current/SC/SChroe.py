################
################
# Schroe.py v0.5
################
################
#
# https://github.com/heedmane/schroepy/
# Licensed under GPLv2 (http://www.gnu.org/licenses/gpl-2.0.html)
#
# Python port of the Mathematica script specified in arXiv:hep-ph/9811453
#
# Author:
# Hector E. Martinez, 
# Physik-Department T30f,
# TU Muenchen,
# Garhing, Germany.
# hector.martinez@tum.de
#
# cfunctions.c include code written by Thomas Rosenhammer
# SClib licensed under GPL v2 https://github.com/drestebon/sclib/



import numpy as np
import matplotlib.pyplot as plt
import scipy as sp
import scipy.integrate as integrate
import math
from scipy.interpolate import pchip


#######
# SClib
#######
# to evaluate C functions within python

import ctypes as ctypes
import sys
import os


##################################################################################################################################
##################################################################################################################################

def error(message):
	sys.stderr.write("error: %s\n" % message)
	sys.exit(1)

class SClib:
	def __init__(self, lib_path, fnames=[]):
        	self.lib_path = lib_path

		self.lib = ctypes.CDLL(lib_path)
		self.TYPE = [ctypes.c_int, ctypes.c_long, ctypes.c_float, ctypes.c_double, ctypes.c_longdouble]
		self.N_INPUTS = {}
		self.INPUT_LEN = {}
		self.INPUT_TYPE = {}
		self.INPUT_DTYPE = {}
		self.N_OUTPUTS = {}
        	self.OUTPUT_LEN = {}
        	self.OUTPUT_TYPE = {}
        	self.OUTPUT_DTYPE = {}

        	self.fnames = []

		for x in fnames:
			if hasattr(self.lib, x):
				self.N_INPUTS[x] = ctypes.c_int.in_dll(self.lib,"_"+x+"_N_INPUTS_").value
				foo = (ctypes.c_int*self.N_INPUTS[x]).in_dll(self.lib,"_"+x+"_INPUT_LEN_")
				self.INPUT_LEN[x] = [foo[i] for i in range(self.N_INPUTS[x])]
				foo = (ctypes.c_int*self.N_INPUTS[x]).in_dll(self.lib,"_"+x+"_INPUT_TYPE_")
				self.INPUT_TYPE[x] = [(self.TYPE[foo[i]]*self.INPUT_LEN[x][i]) for i in range(self.N_INPUTS[x])]
				self.INPUT_DTYPE[x] = [self.TYPE[foo[i]] for i in range(self.N_INPUTS[x])]

				self.N_OUTPUTS[x] = ctypes.c_int.in_dll(self.lib,"_"+x+"_N_OUTPUTS_").value
				foo = (ctypes.c_int*self.N_OUTPUTS[x]).in_dll(self.lib,"_"+x+"_OUTPUT_LEN_")
				self.OUTPUT_LEN[x] = [foo[i] for i in range(self.N_OUTPUTS[x])]
				foo = (ctypes.c_int*self.N_OUTPUTS[x]).in_dll(self.lib,"_"+x+"_OUTPUT_TYPE_")
				self.OUTPUT_TYPE[x] = [(self.TYPE[foo[i]]*self.OUTPUT_LEN[x][i]) for i in range(self.N_OUTPUTS[x])]
				self.OUTPUT_DTYPE[x] = [self.TYPE[foo[i]] for i in range(self.N_OUTPUTS[x])]
				self.fnames.append(x)
			else:
				print("__init__(): "+x+"() discarded: no symbol found")
	def retype(self):
		for x in self.fnames:
			foo = (ctypes.c_int*self.N_INPUTS[x]).in_dll(self.lib,"_"+x+"_INPUT_TYPE_")
			self.INPUT_TYPE[x] = [(self.TYPE[foo[i]]*self.INPUT_LEN[x][i]) for i in range(self.N_INPUTS[x])]
			foo = (ctypes.c_int*self.N_OUTPUTS[x]).in_dll(self.lib,"_"+x+"_OUTPUT_TYPE_")
			self.OUTPUT_TYPE[x] = [(self.TYPE[foo[i]]*self.OUTPUT_LEN[x][i]) for i in range(self.N_OUTPUTS[x])]
	
	def reload(self):
		self.unload()
		self.__init__(self.lib_path,self.fnames)
	
	def unload(self):
		while self.isLoaded():
			try:
				handle = self.lib._handle
			except:
				error("unload(): no se pudo del self.lib")
			libdl = ctypes.CDLL("libdl.so")
			libdl.dlclose(handle)
	
	
	def isLoaded(self):
		libp = os.path.abspath(self.lib_path)
		ret = os.system("lsof -p %d | grep %s > /dev/null" % (os.getpid(), libp))
		return (ret == 0)
	
	def eval(self,fun,*args):
		if not hasattr(self.lib, fun) or fun not in self.fnames:
			error("eval(): "+fun+"() not available")
		if(len(args)!=self.N_INPUTS[fun]):
			error("feval("+fun+"): len(args) = {0} != N_INPUTS = {1}".format(len(args),self.N_INPUTS[fun]))
		py = [self.OUTPUT_TYPE[fun][i]() for i in range(self.N_OUTPUTS[fun])]
		px = [self.INPUT_TYPE[fun][i](*args[i]) for i in range(self.N_INPUTS[fun])]

        	fargs = py+px
		getattr(self.lib, fun)(*fargs)
		return [np.frombuffer(py[i],dtype=self.OUTPUT_DTYPE[fun][i]) for i in range(self.N_OUTPUTS[fun])]


##################################################################################################################################
##################################################################################################################################


#############
# Definitions 
#############


h = 1.e-4   #integration step
ERES = 1e-15 #defines the precision in the eigenvalues

# Definition of the potential function, for instance the Cornell potential, with fixed
# parameters if you want make the potential depend on extra parameters the other functions
# must change accordingly, see example1.

AA = 0.1
TT = 0.5
MM = 1

#the potential
def V(r):
    return (AA/r)+TT*r

#the effective potential
def Veff(r,l):
    return MM*V(r) + l*(l+1)/(r*r)



#####################
# Auxiliary Functions 
#####################

# the following function search for the minimun of the effective potential most to the right 
# always include the point when working with float numbers: 2 is not the same as 2.

def min_(l):
    # looks for the min most to the right (see hep-ph/9811453)
    xrat = 2.   #initial x 
    weit = 0.01 #step
    de  = h/10. #initial step

    aux=[]
    xs = xrat
    if xs < de:
            print 'xmin too small'
    if xs < (de+weit):
	    return de
    rs = Veff(xs+weit,l)	
    ms = Veff(xs,l)	
    ls = Veff(xs-weit,l)	
    while weit > h:
    	
	if xs < (de + weit):
	    return de
	elif rs > ms and ms > ls:
	    xs = xs-weit
	    rs = ms
	    ms = ls
	    ls = Veff(xs-weit,l)
	
	elif rs<ms and ms < ls:
	    xs = xs+weit 
	    ls = ms
	    ms = rs
	    rs = Veff(xs+weit, l)
	elif ls == rs:
	    return xs
	
	elif rs > ms and ls > ms:
	    weit = weit*0.1			
	    rs   = Veff(xs+weit, l)
	    ms   = Veff(xs,l)
	    ls   = Veff(xs-weit,l)
    return xs



# normalize wavefunctions
def normalize(x,y):
    y = np.array(y)
    x = np.array(x)
    y = 1./(math.sqrt(integrate.simps(y*y,dx=h)))*y
    ynorm = np.dstack((x,y))
    return ynorm[0]

# used to an specific range in wfplot2
def simm(x, array):
    tol = h*100
    answer = []
    for y in xrange(len(array)):
        if math.fabs(x-array[y]) <= tol:
            answer.append(y)
    return answer


# load C functions
if __name__ == "__main__":
    from subprocess import check_output
	
    libname = "cfunctions.so"
    libpath = os.getcwd()+"/"+libname
    functions = ['wrap','recover_function','calculate_to_recover']
    try:
	enc_lib = SClib(libpath, functions)
    except:
        pass
	
    if str(check_output(["make", libname])).find("is up to date")<0:
	try:
	    enc_lib.reload()
	except:
	    enc_lib = SClib(libpath, functions)


##################################################################################################################################
##################################################################################################################################


print 'SChroe.py v1.5'
print 'The potential is defined as V(r)  '
print 'the effective potential as Veff(r,l).'
print 'the integration step h = ' +repr(h)
print 'the precision in the eigenvalues is '+str(ERES)
print ' '



#################
## CORE FUNCTIONS
#################


def eigenvalueC(elow1, eup1, n, l):
    # fast finding of the eigenvalue using a C function
    
    elow1 = np.array([elow1])
    eup1 = np.array([eup1])
    n = np.array([n])
    l = np.array([l])
    eig = np.longdouble(enc_lib.eval('wrap',elow1, eup1, n, l)[0][0])
    print 'E = '+str(eig)
    return eig


def wavefunctionC(e,l):
    # fast calculation of the wavefunction using a C function
    e = np.array([e])
    l = np.array([l])
    a = np.array([AA])
    t = np.array([TT])
    m = np.array([MM])
    NNYY = enc_lib.eval('calculate_to_recover',e, l,a,t,m)
    enc_lib.OUTPUT_LEN['recover_function'] = NNYY[0]
    enc_lib.retype()
    YY = enc_lib.eval('recover_function',NNYY[0])
    XX = np.arange(h*0.1,NNYY[0]*h,h)
    return np.dstack((XX,YY[0]))[0]
    


def eigenvaluePy(elow1, eup1, n, l):
    # the following function finds the eigenvalue using Python (no C functions)
	
    de = h*0.1
    feh = MM*ERES
    eup = MM*eup1
    elow = MM*elow1

    xmin=min_(l)+h
    ewidth = eup - elow
	
    while ewidth > feh:

	ewidth = eup - elow
	emid = (eup+elow)*0.5
		
        x  = de
	y  = pow(x,l+1)
	yp = (l+1)*pow(x,l)

	yold = 1.
	n0x = 0

	while True:
	    # Runge-Kutta 4

	    a1 = h*yp
	    b1 = h*(Veff(x,l)-emid)*y

	    a2 = h*(yp+b1*0.5)
	    Vxh = Veff(x+h*0.5,l)-emid
	    b2 = h*Vxh*(y+a1*0.5)

	    a3 = h*(yp+b2*0.5)
	    b3 = h*Vxh*(y+a2*0.5)

	    x = x+h
	    a4 = h*(yp+b3)
	    Vxnew = (Veff(x,l)-emid)
	    b4 = h*Vxnew*(y+a3)

	    y = y + a1/6. + a2/3. + a3/3. +a4/6.
	    yp = yp + b1/6. + b2/3. + b3/3. + b4/6.
			
	    if y*yold < 0:
		n0x=n0x+1

	    yold = y

	    if n0x > n:
		eup = emid # enough nodes
		break

	    if Vxnew > 0 and x > xmin and y*yp > 0: 
		# Veff positive and x > last min and y divergent
		elow = emid
		break
    eig = emid/MM

    print 'E = '+str(eig)
    return eig




def wavefunctionPy(e,l):
    # given the eigenvalue, the following function computes the reduced wave function and return it as a NumPy array using only Python code
	
    xmin=np.longdouble(min_(l)+h)

    ep = np.longdouble(MM*e)
	
    y = []
    r = []

    # beginning of integration

    x  = np.longdouble(h*0.1)
    y.append(np.longdouble(pow(x,l+1)))
    yp = np.longdouble((l+1)*pow(x,l))
    r.append(np.longdouble(h*0.1))

    i=0

    while True:
	# Runge-Kutta of order 4 

	a1 = np.longdouble(h*yp)
	b1 = np.longdouble(h*(Veff(x,l)-ep)*y[i])

	a2 = np.longdouble(h*(yp+b1*0.5))
	Vxh = np.longdouble(Veff(x+h*0.5,l)-ep)
	b2 = np.longdouble(h*Vxh*(y[i]+a1*0.5))

	a3 = np.longdouble(h*(yp+b2*0.5))
	b3 = np.longdouble(h*Vxh*(y[i]+a2*0.5))

	x = np.longdouble(x+h)
	a4 = np.longdouble(h*(yp+b3))
	Vxnew = np.longdouble((Veff(x,l)-ep))
	b4 = np.longdouble(h*Vxnew*(y[i]+a3))

	y.append(np.longdouble(y[i] + a1/6. + a2/3. + a3/3. + a4/6.))
	yp = np.longdouble(yp + b1/6. + b2/3. + b3/3. + b4/6.)
	r.append(x)
	i=i+1

	if Vxnew > 0. and x > xmin and np.longdouble(y[i]*yp) > 0.: 
	    # Veff positive and x > last min and y divergent
	    break
		
    u = normalize(r,y)
    return u




def solve_schroedingerPy(elow,eup,n,l):
    #Finds the eigenvalue and then computes the wavefunction
    eigg = np.array([eigenvaluePy(elow,eup,n,l)])
    yyy = wavefunctionPy(eigg[0],l)
    return eigg,yyy



def solve_schroedingerC(elow,eup,n,l):
    #Finds the eigenvalue and then computes the wavefunction
    eigg = np.array([eigenvalueC(elow,eup,n,l)])
    yyy = wavefunctionC(eigg[0],l)
    return eigg,yyy




###################
## EXTRA FUNCTIONS
###################

# These functions provide like plotting and integration. 
# They don't have a very sofisticated implementation but they could be useful
# at least to perform some basic analysis of the wavefunctions.


def wfplot(ynL):
    #plots a single array (wavefunction) of the form [x,y(x)] in the whole range where it is defined.
    return plt.plot(ynL[:,0],ynL[:,1])



def wfplot2(wflist,bounds):
    # improved ploting, it can plot more than one wavefunction, also it can change the plot range, 
    # if the bounds parameter is -1 it plots the wavefunctions in their full range. It also shows the axis x line.

    # usage: wfplot2([ynL],[0,10]) (NOTE THE BRAKETS) to plot ynL in the interval (0,10) or wfplot2([ynL],-1) to plot ynL in its full range, or 
    # wfplot2([ynL,ynL2],-1) to plot two functions in te full range.

    if bounds == -1:
        plt.axhline(0,color='black')
        for element in wflist:
            wfplot(element)
        return plt.show()
                
    elif len(bounds) == 2:
        plt.axhline(0,color='black')
        for element in wflist:
            xmin = min(simm(bounds[0],element[:,0]))
            xmax = max(simm(bounds[1],element[:,0]))
            xaux = element[:,0][xmin:xmax]
            yaux = element[:,1][xmin:xmax]
            plt.plot(np.array(xaux),np.array(yaux)) 
        return plt.show()

    else:
        print 'invalid bounds'


def wfpro(function, ynL, ynpLp):
    # returns, as an array, the product of two wavefunctions and a continuos function as an array [x,ynL(x)*f(x)*ynpLp(x)] 
    # from the common min value of x to the point corresponding to the min between the two larger values of x.
    aux = []
	
    if len(ynL[:,0]) > len(ynpLp[:,0]):
	ynL = ynL[:,1][:len(ynpLp[:,0])]
	ynL = ynL*ynpLp[:,1]
	for element in ynpLp[:,0]:
	    element = function(element)
	    aux.append(element)
	ynL = ynL*aux
	ynL = zip(ynpLp[:,0],ynL)
	ynL = np.array(ynL)
	return ynL

    elif len(ynL[:,0]) < len(ynpLp[:,0]):
	ynpLp = ynpLp[:,1][:len(ynL[:,0])]
	ynpLp = ynpLp*ynL[:,1]
	for element in ynL[:,0]:
	    element = function(element)
	    aux.append(element)
	ynpLp = ynpLp*aux
	ynpLp = zip(ynL[:,0],ynpLp)
	ynpLp = np.array(ynpLp)
	return ynpLp

    elif len(ynL[:,0]) == len(ynpLp[:,0]):
	ynpLp = ynpLp[:,1]*ynL[:,1]
	for element in ynL[:,0]:
	    element = function(element)
	    aux.append(element)
	ynpLp = ynpLp*aux
	ynpLp = zip(ynL[:,0],ynpLp)
	ynpLp = np.array(ynpLp)
	return ynpLp


def wfpro1(fun, ynLp):
    #returns the product of an array (maybe the product of two functions) and a continues function as an array [x,ynL(x)*f(x)*ynpLp(x)]
    ynL = ynLp.copy()
    aux = []
    for element in ynL[:,0]:
	aux.append(fun(element))
    
    aux = np.array(aux)
    ynL[:,1] = ynL[:,1]*ynL[:,1]
    ynL[:,1] = aux*ynL[:,1]
    return np.array(ynL)
	
	
def wfpro2(ynL, ynpLp):
    #returns the product of two arrays as an array [x,ynL(x)*ynpLp(x)]
    if len(ynL[:,0]) > len(ynpLp[:,0]):
	ynL = ynL[:,1][:len(ynpLp[:,0])]
	ynL = ynL*ynpLp[:,1]
	ynL = zip(ynpLp[:,0],ynL)
	ynL = np.array(ynL)
	return ynL

    elif len(ynL[:,0]) < len(ynpLp[:,0]):
	ynpLp = ynpLp[:,1][:len(ynL[:,0])]
	ynpLp = ynpLp*ynL[:,1]
	ynpLp = zip(ynL[:,0],ynpLp)
	ynpLp = np.array(ynpLp)
	return ynpLp
	
    elif len(ynL[:,0]) == len(ynpLp[:,0]):
	ynpLp = ynpLp[:,1]*ynL[:,1]
	ynpLp = zip(ynL[:,0],ynpLp)
	ynpLp = np.array(ynpLp)
	return ynpLp


def wfpro2a(ynLp):
    #returns the array [x,ynL(x)^2]
    ynl = ynLp.copy()
    ynl[:,1] = ynl[:,1]*ynl[:,1]
    return ynl

def interpolarize(ynL):
    #converts an array of the type [x,y(x)] into a continuous interpolated function f(x) between the points xmin and xmax of the array
    print 'the min value of x is ' + repr(ynL[0][0])
    print 'the max value of x is ' + repr(ynL[-1][0])
    return pchip(ynL[:,0],ynL[:,1])


def wfint(yyf):
    #Simpsons rule integral for an array [x,y(x)]
    return integrate.simps(yyf[:,1],dx=h)


def gaussint(f,a,b,n):
    #Gauss-Legendre quadrature for a continious, maybe interpolated, function
    xw=sp.special.orthogonal.p_roots(n)
    aux=[]
    for x in xw[0]:
	aux.append(((b-a)*0.5)*f(((b-a)*0.5)*x+(a+b)*0.5))
	
    aux = np.array(aux)
    aux=aux*xw[1]
    np.array(aux)
    return np.sum(aux)

def der(R):
    # the following function computes the derivative of an array-defined function [x,y(x)] returning [x,y'(x)]
    Rp = np.gradient(R[:,1],h)
    Rp = zip(R[:,0],Rp)
    Rp = np.array(Rp)
    return Rp


def Rr(y):
    # this function transform the reduced wave radial wave function into the radial wave funcion ( R(r)=y(r)/r )
    R = y[:,1]/y[:,0]
    R = zip(y[:,0],R)
    R = np.array(R)
    return R

