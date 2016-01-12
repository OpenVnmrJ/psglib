
/* tncosyps - homonuclear correlation   (phase-sensitive version)

 Parameters:

      p1 = 90 degree pulse on the observe nucleus
      pw = x degree pulse on the observe nucleus 
           x:  90 degrees for standard COSY experiment

  satdly = length of saturation period;
  sspul  = 'y' for trim(x)-trim(y) before d1
      ss = number of steady state pulses; if ss < 0 then -ss
           steady-state transients are performed before EVERY increment;
           if ss > 0, then ss steady-state transients are performed
           before only the first increment.

   phase =   0:  non-phase-sensitive experiment  (P-type peak selection)
           1,2:  phase-sensitive HYPERCOMPLEX experiment
             3:  phase-sensitive TPPI experiment
      nt = (min):  multiple of 4   (phase = 0)
                   multiple of 2   (phase = 1,2 or 3)
           (max):  multiple of 8   (phase = 0)
                   multiple 0f 4   (phase = 1,2 or 3)


 (all phase cycling taken from s.l. patt, to be submitted)

 s.l.patt   15 may  1985
 revised    24 february  1988  (s.f.)
             2 july      1992   (g.g) */


#include <standard.h>

pulsesequence()
{
/* VARIABLE DECLARATION */
   double          satdly,
                   ss,
                   satpwr,
                   phase;
   int             iphase;
   char            sspul[MAXSTR],satmode[MAXSTR];


/* INITIALIZE VARIABLES */
   ss = getval("ss");
   satdly = getval("satdly");
   phase = getval("phase");
   iphase = (int) (phase + 0.5);
   getstr("satmode",satmode);
   getstr("sspul",sspul);
   satpwr = getval("satpwr");
   if (iphase == 3)
      initval((double) (ix - 1), v14);
   if (p1 == 0.0)
      p1 = pw;


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
      hlv(ct, v9);
      mod4(ct, v2);
      mod2(ct, v1);
   elsenz(ssctr);
      sub(ssval, ssctr, v12);	/* v12 = 0,...,ss-1 */
      hlv(v12, v9);
      mod4(v12, v2);
      mod2(v12, v1);
   endif(ssctr);


/* CALCULATE PHASECYCLE */
   hlv(v9, v8);
   add(v8, v9, v9);
   mod2(v9, v9);		/* 00111100 */
   dbl(v1, v1);
   add(v1, v9, v1);		/* 0202+00111100 */
   initval(4.0, v10);
   sub(v10, v2, v2);		/* 0321 */
   if ((iphase == 1) || (iphase == 2))
      assign(zero, v2);
   add(v2, v9, v2);		/* 0321+00111100 or 0+00111100 */
   mod4(v1, oph);
   if (iphase == 2)
      incr(v1);
   if (iphase == 3)
      add(v1, v14, v1);


/* BEGIN ACTUAL PULSE SEQUENCE */
   status(A);
     if (sspul[A]=='y')
      {
       rlpower(tpwr-12,TODEV);
       rgpulse(200*pw,zero,rof1,0.0);
       rgpulse(300*pw,one,0.0,rof1);
       rlpower(tpwr,TODEV);
      }
     if (satmode[A]=='y')
      {
       rlpower(satpwr,TODEV);
       rgpulse(satdly,zero,rof1,rof2);
       rlpower(tpwr,TODEV);
      }
     hsdelay(d1);
   status(B);
      rgpulse(p1, v1, rof1, 1.0e-6);
     if (satmode[B]=='y')
      {
      if (d2>0)
       {
        rlpower(satpwr,TODEV);
        rgpulse(d2-3.0*rof1-9.4e-6-4*pw/3.1414,zero,rof1,rof1);
        rlpower(tpwr,TODEV);
       }
      }
     else
      delay(d2 -rof1 -1.0e-6 -4*pw/3.1414);
      rgpulse(pw, v2, rof1, rof2);
   status(C);
}
