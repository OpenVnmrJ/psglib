/* VARIAN VXR-5000 MENU GENERATED PULSE SEQUENCE:  hahnse */

/* HAHN SPIN ECHO: Hahn, E.L.,Phys. Rev. 80,580-594,1950. */

/* ....D1....[pi/2=PW]....D2....[2PW]....D2....ACQ */

/* ....A.....[...................B............][.C.] */

/*Fernando Commodari, January 15,1989; revised April 29, 1992*/
 
#include <standard.h>
 
pulsesequence()
{
initval(1.0,v9); divn(ct,v9,v8); /* v8=Cyclops */
status(A); 
delay(d1);
status(B); 
add(zero,v8,v1); pulse(pw,v1);
delay(d2);
add(zero,v8,v1); pulse(2.0*pw,v1);
delay(d2);
status(C); 
}
