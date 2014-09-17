/*################
################
# Schroe.py v0.5
################
################
#
# https://github.com/heedmane/schroepy/
# Licensed under GPL v2 (http://www.gnu.org/licenses/gpl-2.0.html)
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
*/




#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <sclib.h>
#include "potential.h"
#include <string.h>

long double * __RECOVERED_WAVEFUNCTION__;

typedef struct{
	// this structures stores a wavefunction array that nows its number of elements N
	// this is necessary to use the integration functions with no major changes
	int N; //size of the array
	long double *y; //wavefunction
} wf;



long double minimum(int l){
	long double xrat=2, weit=0.01, del=h/10;
	long double rs, ls, ms, xs=xrat;
	
	
	if(xs<del){
		printf("xmin too small\n");
		return 0;
	}
	
	if(xs<del+weit){return del;}
		
		rs=Veff(xs+weit,l);
		ms=Veff(xs,l);
		ls=Veff(xs-weit, l);
	
	while(weit>h){
		if(xs<del+weit){return del;}
		
		if((rs>ms)&&(ms>ls)){
			xs=xs-weit;
			rs=ms;
			ms=ls;
			ls=Veff(xs-weit,l);			
			}
		else if((rs<ms)&&(ms<ls)){
			xs=xs+weit;
			ls=ms;
			ms=rs;
			rs=Veff(xs+weit,l);
			}
		else if(ls==rs){
			return xs;}
		else if((rs>ms)&&(ls>ms)){
			weit=weit/10;
			rs=Veff(xs+weit,l);
			ms=Veff(xs, l);
			ls=Veff(xs-weit,l);
		}
		
	}
	
	return xs;
}

long double Veffminusepsilon(long double r, int l, long double epsilon){
	return Veff(r,l)-epsilon;
}


long double eigenvalue(long double elow1, long double eup1, int n, int l){
	long double del=h/10, feh=ERES*MM, eup=MM/hbar/hbar*eup1, elow=MM/hbar/hbar*elow1, xmin, x;
	long double ewidth, emid, y, yp, yold;
	long double a1, a2, a3, a4, b1, b2, b3, b4, Vxh, Vxneu=-1;
	int n0x;

		
	xmin=minimum(l)+h;
	
	do{
		ewidth=eup-elow;
		emid=(eup+elow)/2;
		
		//Anfang Integration
		x=del;
	
			y=powl(x,l+1);
			yp=(l+1)*powl(x,l);
		
		yold=1;
		n0x=0; //Anzahl der Knoten
		do{
		
			// nächstes y Runge-Kutta
		
			a1=h*yp;
			b1=h*Veffminusepsilon(x,l, emid)*y;
		
		
			a2=h*(yp+b1/2);
			Vxh=Veffminusepsilon(x+h/2, l,emid);
			b2=h*Vxh*(y+a1/2);
		
			a3=h*(yp+b2/2);
			b3=h*Vxh*(y+a2/2);
		
			x=x+h;
			a4=h*(yp+b3);
			Vxneu=Veffminusepsilon(x,l, emid);
			b4=h*Vxneu*(y+a3);
		
			y=y+a1/6+ a2/3+a3/3+a4/6;
			yp=yp+b1/6+ b2/3+b3/3+b4/6;
		
			// Auf Knoten prüfen:
			if(y*yold<0){
				n0x++;} //neuer Knoten
				
			yold=y;
				
			if(n0x>n){
				eup=emid; // zuviele Knoten
				break;}
			if(((Vxneu>0)&&(x>xmin))&&(y*yp>0)){ // Veff positiv + x > letztes Minimum + y divergiert
				elow=emid;
				break;}
		}while(1); // immer erfüllt
		
	}while( ewidth>feh);
	
	return emid/MM*hbar*hbar;
}





long double integration(int N, long double f[N]){
	int i;
	long double Th, T2h, T4h, abl=f[N-1]-f[N-2]-f[0], abl3=f[N-1]-3*f[N-2]+3*f[N-3]-f[N-4]-3*f[0]+3*f[1]-f[2];
	
	Th=f[N-1]/2; // f an der Stelle 0 = 0 da u = r^(l+1) für r->0;
	T2h=Th;
	T4h=Th;
	
	for(i=0;i<=N-2;i++){
		Th+= f[i];
		if(2*i+1<=N-3){
			T2h+=f[2*i+1];
		}
		if(4*i+3<=N-5){
			T4h+=f[4*i+3];
		}
	
	}
	
	Th=h*(Th-abl/12.-abl3/24./30.);
	T2h=2*h*(T2h-abl/6.-abl3/90.);
	T4h=4*h*(T4h-abl/3.-abl3*4/45.);


	return (64*Th-20*T2h+T4h)/45.;
	
}




wf wavefunction(long double e,int L,long double a,long double t,long double m){

	wf result; // to store the result

	// variables to use in the RK
	long double xmin, ep, x, yp;
	long double a1, b1, a2, b2, a3, b3, a4, b4, Vxh, Vxnew;

	result.y = malloc(sizeof(long double)*5);
        if (result.y==NULL){
                printf("Error allocating memory!\n"); 
                return result; 
        }

	xmin = minimum(L);

	ep = m*e;
	x = h;

	//begin of the integration according to initial conditions
	result.y[0] = powl(x,L+1);
	yp = (L+1)*powl(x,L);

	
	int i = 0;

	// Runge-Kutta order 4
	do{
		
		a1 = h*yp;
		b1 = h*Veffminusepsilon(x,L,ep)*result.y[i];

		a2 = h*(yp+b1*0.5);
		Vxh = Veffminusepsilon(x+h*0.5,L,ep);
		b2 = h*Vxh*(result.y[i]+a1*0.5);

		a3 = h*(yp+b2*0.5);
		b3 = h*Vxh*(result.y[i]+a2*0.5);
		a4 = h*(yp+b3);
		Vxnew =  Veffminusepsilon(x+h,L,ep);
		b4 = h*Vxnew*(result.y[i]+a3);
	


                long double *temp = realloc(result.y, (i+2)*sizeof(long double));
                if ( temp != NULL ) //realloc was successful
                {
                    result.y = temp;
                    //printf("more memory given!\n");
                }
                else //there was an error
                {
                    free(result.y);
                    printf("Error allocating memory!\n");
                    return result;
                }


		result.y[i+1] = result.y[i] + a1/6. + a2/3. + a3/3. + a4/6.;
		yp = yp + b1/6. + b2/3. + b3/3. + b4/6.;

		//printf("y = %Lf\n",result.y[i]);
		
		x = x+h;
		i = i+1;

	}while(!(Vxnew > 0 && x > xmin && result.y[i]*yp > 0 ));
	
	result.N = i;
	int j;

	long double u2[result.N];

	for(j=0;j<result.N;j++){
		u2[j] = result.y[j]*result.y[j];
	}	

	long double norm = integration(result.N,u2);
	
	for(j=0;j<result.N;j++){
		result.y[j] = (1/sqrtl(norm))*result.y[j];
	}

	return result;

}

/* functions that will be called in the python script */

PYO(wrap,1, 1);
PYO_TYPES(wrap,1, LDOUBLE);
PYI(wrap,4,1,1, 1, 1);
PYI_TYPES(wrap,7, LDOUBLE, LDOUBLE, INT, INT);
void wrap(long double * result, long double * elow, long double * eup, int * n, int * l){

    *result =  eigenvalue(elow[0],eup[0],n[0],l[0]);

}



PYO(recover_function,1, 0);
PYO_TYPES(recover_function,1, LDOUBLE);
PYI(recover_function,1,1);
PYI_TYPES(recover_function,1, INT);
void recover_function(long double * OO, int *SIZE){
    // recovers the wavefunction from the functions that looks for enigenvalue
    memcpy(OO, __RECOVERED_WAVEFUNCTION__, SIZE[0]*sizeof(long double));

}


PYO(calculate_to_recover,1, 1);
PYO_TYPES(calculate_to_recover,1, INT);
PYI(calculate_to_recover,5,1,1,1,1,1);
PYI_TYPES(calculate_to_recover,5, LDOUBLE, INT,LDOUBLE,LDOUBLE,LDOUBLE);
void calculate_to_recover(int *result,long double * e,int *L,long double *a,long double *t,long double *m){

    wf NNYY;
    NNYY =  wavefunction(e[0],L[0],a[0],t[0],m[0]);
    __RECOVERED_WAVEFUNCTION__ = NNYY.y;
    result[0] = NNYY.N;

}



