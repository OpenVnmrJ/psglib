/* VARIAN VXR-5000  PULSE SEQUENCE:  tqfb */

/* TO DETECT THE TRIPLE QUANTUM COHERENCE ARISING FROM THE */
/* BIEXPONENTIAL TRANSVERSE RELAXATION OF SODIUM */
/* Chun-wa and Stephen Wimperis, J. Magn. Res. 88 440,1990*/
/* Personal Communication, April 30, 1990*/
/*Fernando Commodari, May 29 1990; revised April 29, 1992.*/

/*
- d1 - p1 (t1) - d2 - p2 (t2) - d2 - p1 (t3) - d3 - p1 (t4) - acq (t5)

	d1 = optional relaxation delay

	p1 = pi/2 pulse with phase t1(echo) t3 (TQ creation) and t4 (read pulse)

	d2 = 1/2 the echo (TQ preparation) time (msec); 
            not < (rof1+rof2+p2/2+p1/2)

	p2 = 180 pulse with phase t2

	d3 = TQ evolution delay (microsec); not < (rof2+rof1+p1)

	t5 = receiver phase */

#include <standard.h>
 
pulsesequence()
{
double d3, p2;
d3=getval("d3");
p2=getval("p2");

/*PERIOD 1 */;
 status(A); 
 delay(d1);
 loadtable("tqfb");
 stepsize(30.0,TODEV); 

/*PERIOD 2 */

 status(B);
 xmtrphase(t1);
 pulse(p1,zero);
 delay(d2-p1/2-p2/2-rof1-rof2);
 xmtrphase(t2);
 pulse(p2,zero);
 delay(d2-rof1-rof2-p2/2-p1/2);
 xmtrphase(t3); 
 pulse(p1,zero);
 delay(d3-rof2-rof1-p1);

/*PERIOD 3*/
 status(C);
 stepsize(90.0,TODEV);
 xmtrphase(t4);
 pulse(p1,zero);
 setreceiver(t5);   

 }
