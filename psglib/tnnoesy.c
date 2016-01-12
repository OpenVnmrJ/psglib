
/* tnnoesy -  2D cross relaxation experiment.  It can be performed in
            either phase-sensitive or absolute value mode.  Either
            TPPI or the hypercomplex method can be used to achieve
            F1 quadrature in a phase-sensitive presentation.  No
            attempt is made to suppress J-cross peaks in this pulse
            sequence.
                TRANSMITTER SOLVENT SATURATION ONLY
                ASSUMES ON-RESONANCE SOLVENT (tof is at solvent position)

   satmode  determines when the saturation happens. Satmode should be set
            analogously to dm, i.e. satmode='yyyn' or 'ynyn' etc.
               satmode='yyyn' recommended
    satdly = length of presaturation period (saturation may also occur in d2 and mix
              as determined by satmode)
        ss = number of steady state pulses; if ss < 0 then -ss
             steady-state transients are performed before EVERY increment;
             if ss > 0, then ss steady-state transients are performed
             before only the first increment.
     sspul = 'y': selects for Trim(x)-Trim(y) sequence at start of pulse sequence
              G. Gray Palo Alto  Sept. 1991
*/


#include <standard.h>

pulsesequence()
{
   double          arraydim,
                   ss,
                   satdly,satpwr,
                   ni,
                   mix;
   int             iphase;
   char            sspul[MAXSTR],satmode[MAXSTR];


/* LOAD VARIABLES */
   ni = getval("ni");
   arraydim = getval("arraydim");
   satdly = getval("satdly");
   satpwr = getval("satpwr");
   getstr("satmode",satmode);
   mix = getval("mix");
   iphase = (int) (getval("phase") + 0.5);
   ss = getval("ss");
   getstr("sspul", sspul);
   initval(satpwr,v11); initval(tpwr,v13);

   if (iphase == 3)
      initval((double)((int)((ix-1)/(arraydim/ni)+1.0e-6)), v14);
   else
      assign(zero, v14);


/* CHECK CONDITIONS */
   if ((rof1 < 9.9e-6) && (ix == 1))
      fprintf(stdout,"Warning:  ROF1 is less than 10 us\n");

   if (satpwr > 40)
        {
         printf("satpwr too large  - acquisition aborted./n");
         abort(1);
        }

/* DETERMINE STEADY-STATE MODE */
   if (ss < 0)
   {
      ss = (-1) * ss;
   }
   else
   {
      if ((ss > 0) && (ix == 1))
      {
	 ss = ss;
      }
      else
      {
	 ss = 0;
      }
   }
   initval(ss, ssctr);
   initval(ss, ssval);


/* STEADY-STATE PHASECYCLING
/* This section determines if the phase calculations trigger off of (SS - SSCTR)
   or off of CT */

   ifzero(ssctr);
      mod2(ct, v2);
      hlv(ct, v3);
   elsenz(ssctr);
      sub(ssval, ssctr, v12);	/* v12 = 0,...,ss-1 */
      mod2(v12, v2);
      hlv(v12, v3);
   endif(ssctr);


/* CALCULATE PHASECYCLE */
   dbl(v2, v2);
   hlv(v3, v10);
   hlv(v10, v10);
   if (iphase == 0)
   {
      assign(v10, v9);
      hlv(v10, v10);
      mod2(v9, v9);
   }
   else
   {
      assign(zero, v9);
   }
   assign(v10,v1);
   hlv(v10, v10);
   mod2(v1, v1);
   dbl(v1, v1);
   add(v9, v2, v2);
   mod2(v10, v10);
   add(v1, v2, oph);
   add(v3, oph, oph);
   add(v10, oph, oph);
   add(v10, v1, v1);
   add(v10, v2, v2);
   add(v10, v3, v3);
   add(v10,v14,v5);
   if (iphase == 2)
      { incr(v2); incr(v5); }
   if (iphase == 3)
      add(v2, v14, v2);		/* TPPI phase increment */

/*HYPERCOMPLEX MODE USES REDFIELD TRICK TO MOVE AXIAL PEAKS TO EDGE */
    initval(2.0*(double)(((int)(d2*getval("sw1")+0.5)%2)),v6);
  if ((iphase==1)||(iphase==2))
       {add(v2,v6,v2); add(oph,v6,oph); add(v5,v6,v5);}  

/* BEGIN THE ACTUAL PULSE SEQUENCE */
   status(A);
      power(zero,DODEV);
      if (sspul[A] == 'y')
       {
        rgpulse(200*pw, zero, rof1, 0.0e-6);
        rgpulse(200*pw, one, 0.0e-6, rof1);
       }
      if (d1>hst) hsdelay(d1);
      if (satmode[A] == 'y')
      {
        power(v11,TODEV);
        rgpulse(satdly,v5,rof1,rof1);
        power(v13,TODEV); 
      }
   status(B);
      rgpulse(pw, v2, rof1, 1.0e-6);
      if (satmode[B] =='y')
       {  
        if (d2 > 0.0) 
         { 
          power(v11,TODEV);
          rgpulse(d2 - 9.4e-6 - rof1 - 10.0e-6 - (4.0*pw/3.14159),zero,5.0e-6,5.0e-6);
          power(v13,TODEV);
         }
       }
      else
       {
        if (d2 > 0.0)
         delay(d2 - 1.0e-6 - rof1 - (4.0*pw/3.14159));
       }
      rgpulse(pw, v1, rof1, 1.0e-6);
   status(C);
      if (satmode[C] == 'y')
       {
          hsdelay(hst);
          power(v11,TODEV);
          rgpulse(mix-hst,zero,2.0e-6,rof1);
          power(v13,TODEV); 
       }
      else
          hsdelay(mix);
   status(D);
      rgpulse(pw, v3, rof1, rof2);

/*  Phase cycle:   ...satdly(v5)...pw(v2)..d2..pw(v1)..mix...pw(v3)..at(oph)
    (for phase=1. for phase = 2 incr(v2) and incr(v5) )
        v2 =[02] for axial peaks
        v1 =[02]16 for axial peaks
        v3 =[0123]2  "4 step phase cycle selection"
       v10 =[01]8   for quad image
       oph = v1+v2+v3+v10

v5: 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1
v2: 0 2 0 2 0 2 0 2 0 2 0 2 0 2 0 2 1 3 1 3 1 3 1 3 1 3 1 3 1 3 1 3
v1: 0 0 0 0 0 0 0 0 2 2 2 2 2 2 2 2 1 1 1 1 1 1 1 1 3 3 3 3 3 3 3 3
v3: 0 0 1 1 2 2 3 3 0 0 1 1 2 2 3 3 1 1 2 2 3 3 0 0 1 1 2 2 3 3 0 0
oph:0 2 1 3 2 0 3 1 2 0 3 1 0 2 1 3 1 3 2 0 3 1 0 2 3 1 0 2 1 3 2 0    */
}
