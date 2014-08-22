/*################
################
# Schroe.py v1.0
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


PYO(wrap,1, 1);
PYO_TYPES(wrap,1, LDOUBLE);
PYI(wrap,4,1,1, 1, 1);
PYI_TYPES(wrap,7, LDOUBLE, LDOUBLE, INT, INT);
void wrap(long double * result, long double * elow, long double * eup, int * n, int * l){

    *result =  eigenvalue(elow[0],eup[0],n[0],l[0]);

}

