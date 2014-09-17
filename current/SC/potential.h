/* example of a potential.h file 
 * the integration step h must be defined here
 * also hbar ( hbar = 1 for natural units)
 * ERES defines the precision in the eigenvalue
 */


#define hbar 1
#define h 1e-4
#define ERES 1e-15


/* Definition of the potential function, for instance the Cornell potential, with fixed
 parameters. If you want make the potential depend on extra parameters the other functions
 must change accordingly, see example1 */

#define AA 0.1
#define TT 0.5
#define MM 1.


long double V(long double r){
	return (AA/r)+TT*r;
}


long double Veff(long double r, int l){

        return MM*V(r)+l*(l+1)/(r*r);
}


