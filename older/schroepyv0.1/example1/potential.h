/* example of a potential.h file 
 * the integration step h must be defined here
 * also hbar ( hbar = 1 for natural units)
 * ERES defines the precision in the eigenvalue
 */
/* In this example the potential depends on more parameters than r and the functions have been changed accordingly */

#define hbar 1
#define h 1e-4
#define ERES 1e-15

long double V(long double r, long double k, long double sig){
	return (k/r)+sig*r;
}


long double Veff(long double r, long double k, long double sig, long double m, int l){

        return m*V(r,k,sig)+l*(l+1)/(r*r);
}


