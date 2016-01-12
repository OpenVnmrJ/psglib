#ifndef LINT
static char     SCCSid[] = "@(#)roesy1d.c 3.4 8/24/90  Copr 1990 P.S.";
#endif

/* roesy1d - one dimensional rotating frame NOE experiment;
             a continuous-wave spin lock is employed;
             compensation for off-resonance effects may be performed.

   ref:      H. Kessler, U. Anders, G. Gemmecker and S. Steuernagel
             J. Magn. Reson. 85, 1-14 (1989). 

  Parameters:
        p1 = 90 degree pulse on protons (power level at "p1lvl")
     p1lvl = power level for the p1 pulse 
      tpwr = power level for the spin lock pulse
       mix = mixing time
        ss = number of steady state pulses
     sspul = 'y': selects for trim(x)-trim(y) sequence at start
    rocomp = 'n': no resonance offset compensation
             'y': resonance offset compensation
    satdly = length of saturation (at power satpwr); if
    selpwr = transmitter power for the shaped pulse
     selpw = 90 degree shape pulse in microseconds
  selshape = shape of the shaped pulse
    selfrq = transmitter frequency for selective excitation
    satpwr = power for presaturation
    satfrq = frequency for presaturation
    satmode = 'y': selects  presaturation
             'n': no presaturation
     purge = 'y': selects compensating purging pulse
             'n': purging pulse is bypassed
        nt = min:  multiple of 8 
             max:  multiple of 16   (recommended)

  P. Sandor (Darmstadt, modified gg. sept 92 */

#include <standard.h>

pulsesequence()
{
   double     
         satdly = getval("satdly"),
         selpwr = getval("selpwr"),
         mix = getval("mix"),
         selpw = getval("selpw"),
         selfrq = getval("selfrq"),
         p1lvl = getval("p1lvl"),
         satpwr = getval("satpwr"),
         satfrq = getval("satfrq");
   int   
         roc_flag;
   char 
         sspul[MAXSTR],
         rocomp[MAXSTR],
         selshape[MAXSTR],
         satmode[MAXSTR],
         purge[MAXSTR];

/* LOAD AND INITIALIZE PARAMETERS */
         getstr("sspul", sspul);
         getstr("rocomp", rocomp);
         getstr("satmode", satmode);
         getstr("selshape", selshape);
         getstr("purge", purge);

         roc_flag=(rocomp[0] == 'y');

/* CHECK CONDITIONS */
         if (rof1 < 10.0e-6) rof1=10.0e-6;

/* SETUP PHASES */
         assign(one,v1);
         mod2(ct,v11);          /* 0101 */
         hlv(ct,v9);            /* 0011 2233 */
         mod2(v9,v9);           /* 0011 */
         add(v9,v11,v2);        /* 0112 */
         mod2(v2,v2); dbl(v2,v2); /* 0220 */
         dbl(v9,v3);            /* 0022 */
         hlv(ct,oph); hlv(oph,oph); /* 0000 1111 2222 3333    */
         add(v1,oph,v1);  /* 1111 2222 3333 0000   shaped pulse  */
         add(v2,oph,v2);  /* 0220 1331 2002 3113   purging pulse */
         add(v3,oph,v3);  /* 0022 1133 2200 3311   spin-lock     */
         add(one,oph,oph); /*1111 2222 3333 3333   OPH*/

/* BEGIN ACTUAL PULSE SEQUENCE */
   status(A);
      if (sspul[0] == 'y')
      {
         rlpower(p1lvl, TODEV); 
         rgpulse(200*p1, zero, rof1, 0.0);
         rgpulse(220*p1, one, 0.0,rof1);
      }
      hsdelay(d1);
      if (satmode[A] == 'y')
       {
         offset(satfrq,DODEV);
         rlpower(satpwr,TODEV);
         rgpulse(satdly,zero,rof1,rof2);
         offset(tof,TODEV);
         delay(5.0e-6);
       }

   status(B);
      if (selpw!=0.0)
       {
         offset(selfrq,TODEV);
         rlpower(selpwr,TODEV);   
         delay(25.0e-6);
         shaped_pulse(selshape,selpw,v1,rof1,rof2);
         offset(tof,TODEV); 
       }
         rlpower(p1lvl,TODEV);
         delay(25.0e-6); 
     if (purge[A] == 'y')
        rgpulse(p1,v2,rof1,10.0e-6);    /* purging pulse */
     if (roc_flag) rgpulse(p1,v3,rof1,rof2); /* for off. res. compensation */
     if (mix > 0.0)
       {
         rlpower(tpwr,TODEV);  /* for spin-lock */
         rgpulse(mix,v3,rof1,0.0);
         if (roc_flag) 
           {
             rlpower(p1lvl,TODEV); 
             rgpulse(p1,v3,rof1,rof2);
           }
        }

   status(C);
}

