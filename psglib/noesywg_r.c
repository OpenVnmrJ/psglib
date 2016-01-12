
/* noesywg_r -  NOESY experiment with water suppression by gradient echo. 
               Either TPPI or the hypercomplex method can be used to achieve
               F1 quadrature in a phase-sensitive presentation.  No
               attempt is made to suppress J-cross peaks in this pulse
               sequence. F1 axial peaks are shifted by States TPPI method.

               ech dast feb.93  gg palo alto jan 95
               added flipback 16 april 95
               added gradients in t1 to eliminated radiation damping effects
                (gzlvl3 can be set very low, gzlvl3=50 suggested)
               (following suggestion of V.Sklenar, JMR, 114, 132(1995) )

               set hst=0 and gt1=mix to enhance water exchange crosspeaks by
               preventing radiation damping during mix. Phase of flipback
               pulse is currently only optimal for the case of radiation
               damping during mix, so set flipback='n' for enhancing
               water exchange crosspeaks. With proper phase cycling of the
               flipback pulse, it could be used. Again, set gzlvl1=50.
                5jan96 GG
*/
#include <standard.h>
pulsesequence()
{
   double          arraydim,
                   ss,
                   p180 = getval("p180"),
                   gzlvl1 = getval("gzlvl1"),
                   gt1 = getval("gt1"),
                   gzlvl0 = getval("gzlvl0"),
                   gt0 = getval("gt0"),
                   gzlvl2 = getval("gzlvl2"),
                   gt2 = getval("gt2"),
                   gzlvl3 = getval("gzlvl3"),
                   grise = getval("grise"),
                   phincr1 = getval("phincr1"),
                   phincr2 = getval("phincr2"),
                   p1lvl = getval("p1lvl"),
                   mix = getval("mix");
   int             t1_counter,
                   iphase;
   char            flipback[MAXSTR];


/* LOAD VARIABLES */
   getstr("flipback",flipback);
   arraydim = getval("arraydim");
   iphase = (int) (getval("phase") + 0.5);
   ss = getval("ss");
   if (phincr1 < 0.0) phincr1=360+phincr1;
   initval(phincr1,v8);
   if (phincr2 < 0.0) phincr2=360+phincr2;
   initval(phincr2,v11);

   if (iphase == 3)
   {
      t1_counter = ((int) (ix - 1)) / (arraydim / ni);
      initval((double) (t1_counter), v14);
   }
   else
      assign(zero, v14);


/* CHECK CONDITIONS */
   if ((rof1 < 9.9e-6) && (ix == 1))
      fprintf(stdout,"Warning:  ROF1 is less than 10 us\n");


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
      zgradpulse(gzlvl0,gt0);
      obsstepsize(45.0);
      initval(7.0,v7);
      xmtrphase(v7);
      delay(d1);
      rcvroff();
 status(B);
      rgpulse(pw, v2, rof1, rof1);
      if (d2 > 0.0)
       {
     zgradpulse(gzlvl3,0.4*d2-SAPS_DELAY/2.0-GRADIENT_DELAY-(2.0*pw/PI));
     delay(0.1*d2-rof1);
     zgradpulse(-gzlvl3,0.4*d2-SAPS_DELAY/2.0-GRADIENT_DELAY-(2.0*pw/PI));
     delay(0.1*d2-rof1);
       }
      xmtrphase(zero);
      rgpulse(pw, v1, rof1, 0.0);
 status(C);
      delay(hst);
      zgradpulse(gzlvl1,gt1);
      add(v3,two,v4);
      if (flipback[A] == 'y')
       {
        obspower(p1lvl);
        obsstepsize(1.0);
        xmtrphase(v8);
        delay(mix-p1 -hst-gt1);
        rgpulse(p1,v4,20.0e-6,0.0);
        xmtrphase(zero);
        obspower(tpwr);
       }
      else
        delay(mix - hst -gt1);
      rgpulse(pw, v3,30.0e-6,rof2);
     if (p180>0.0) 
      {
      delay(tau-rof2-pw/2);
      zgradpulse(gzlvl2,gt2+grise);
      delay(grise-2.0*SAPS_DELAY);
      obspower(p1lvl);
      obsstepsize(1.0);
      xmtrphase(v11);
      rgpulse(p1,v4,rof1,rof1);
      obspower(tpwr);
      xmtrphase(zero);
      rgpulse(p180, v3,rof1,rof1);
      obspower(p1lvl);
      rgpulse(p1,v4,rof1,rof1);
      obspower(tpwr);
      delay(tau); 
      zgradpulse(gzlvl2,gt2+grise);
      delay(grise);
     }
   status(D);

/* phase cycle: .....pw(v2)..d2..pw(v1)..mix...pw(v3)..at(oph)
    (for phase=1 for phase = 2 incr(v2) and incr(v5) )
     
        v1 =[02]16 for axial peaks
        v3 =[0123]2  "4 step phase cycle selection"
       v10 =[01]8   for quad image
       oph = v1+v2+v3+v10

v2: 0 2 0 2 0 2 0 2 0 2 0 2 0 2 0 2 1 3 1 3 1 3 1 3 1 3 1 3 1 3 1 3
v1: 0 0 0 0 0 0 0 0 2 2 2 2 2 2 2 2 1 1 1 1 1 1 1 1 3 3 3 3 3 3 3 3
v3: 0 0 1 1 2 2 3 3 0 0 1 1 2 2 3 3 1 1 2 2 3 3 0 0 1 1 2 2 3 3 0 0
oph:0 2 1 3 2 0 3 1 2 0 3 1 0 2 1 3 1 3 2 0 3 1 0 2 3 1 0 2 1 3 2 0    */
}
