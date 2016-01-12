/*  sg2pul - shaped gradient two-pulse sequence 

    Parameters: 
        gzlvl1 = gradient amplitude (-32768 to +32768)
        gt1 = gradient duration in seconds (0.002)   
 
    Modified by R.A. Byrd, NCI-FCRDC for shaped gradients 
    and to satisfy dps display(when you use dps the sequence
     defaults to zgradpulse, for display purposes)

*/

#include <standard.h>
#include <shapegrad.c>

extern int dps_flag;

pulsesequence()
{
   double gzlvl1,gt1;

   gzlvl1 = getval("gzlvl1");
   gt1 = getval("gt1");

status(A);
   lk_sample(); 
   delay(d1);
   lk_hold();

status(B);
   rgpulse(p1, zero,rof1,rof2);
   if (dps_flag)
   {
    zgradpulse(gzlvl1,gt1);
   }
   else
   {
    shapegrad(gzlvl1,gt1);
   }
   delay(d2);

status(C);
   pulse(pw,oph);
}

/*
*****************************************************************************
    PUT THE FILE GIVEN BELOW IN YOUR PSGLIB
*****************************************************************************

/* 
 * shapegrad.c
 * 
 * written by R.A. Byrd, Macromolecular NMR Section,
 *            ABL-Basic Research Program
 *            NCI-FDCRDC
 *            Frederick, MD 21702
 *            rabyrd@nmrsgi2.ncifcrf.gov
 *            December 1993
 *
 * this file should be used as an include file within
 * pulsesequences.  It permits the user to call a 
 * shaped gradient by passing the arguments of amplitude
 * and duration.  Adjustments are made in the routine
 * to compensate for dac turn on times as observed 
 * experimentally....they are pretty close, such that
 * when a user asks for a 1.0 msec gradient it comes
 * out that way in real time.
 *
 */
/* calculate gradient ramp shapes */
/* create function for shaped gradient */

shapegrad(gval,gtim)
double gval,gtim;
{
   static double grad_up[8]= {0.0,0.1392,0.3001,0.4705,
			      0.6366,0.7842,0.9003,0.9745},
		 grad_dn[8]= {1.0,0.9745,0.9003,0.7842,
			      0.6366,0.4705,0.3001,0.1392};
   double gstep,gampl;
   int jcnt;
/* determine time step for gradient ramp up */
   if ((0.2*gtim) < (7*8.65e-6))
      gstep = 8.65e-6;
   else 
      gstep = ((0.2*gtim)/7.0);

/* perform the ramped gradient */
    for (jcnt=1;jcnt<=7;jcnt++) {
      gampl=(grad_up[jcnt]*gval);
      rgradient('z',gampl);
      delay(gstep - GRADIENT_DELAY);
      }
    if (gstep > 8.65e-6)
     { rgradient('z',gval);
       delay((0.6*gtim)-(2*GRADIENT_DELAY));}
    else
     { rgradient('z',gval);
       delay(gtim - (14*8.65e-6) - GRADIENT_DELAY);}
/*       delay(1.0e-3 - GRADIENT_DELAY);}*/
    for (jcnt=1;jcnt<=7;jcnt++) {
       gampl=(grad_dn[jcnt]*gval);
       rgradient('z',gampl);
       delay(gstep - GRADIENT_DELAY);
     }
   rgradient('z',0.0);
}

*/
