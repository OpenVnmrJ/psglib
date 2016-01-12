/* selhmqc- selective HMQC with hypercomplex and
         TPPI implementations, steady-state transients for either all or
         just the first t1 increment, pulse+receiver phasecycling during
         steady-state transients, and presaturation using the transmitter
         using the waveform generator

   ref: Summers et al., j. amer. chem. soc. 108:4285-4294 (1986)


   Parameters:
      sspul = 'y' does trim(x)-trim(y) where trim is 200*pw (at tpwr) to
              destroy all 1H magnetization (forces steady-state)

     satmode = 'yn':  presaturation during relaxation period (satdly) with xmtr
               'yy':  presaturation during both relaxation period (satdly)
                      and null period (null) with xmtr
     satfrq = presaturation frequency
     satdly = saturation time during the relaxation period
     satpwr = saturation power for all periods of presaturation with xmtr
     selpwr = selective decoupler pulse power
     selpw  = selective decoupler pulse length (in usec) at dof
   selshape = selective decoupler (wfg) shape
        dof = frequency for decoupling during acquisition
         hs = 'yn':  homospoil pulse (hst) during the d1 relaxation delay
              'yy':  homospoil pulse during both d1 and nulling period (null)
         pw = 90 degree xmtr pulse length for protons (the observed nucleus)
     pwxlvl = power level for decoupler pulses
        pwx = 90 degree decoupler pulse length for the heteronucleus
          j = one-bond heteronuclear coupling constant (in Hz)
       jeff = effective j (longer than j to account for precession during selpw)
         dm = 'nnn':  no broadband decoupling of heteronucleus during
                      acquisition; 
              'nny':  broadband heteronuclear decoupling during acquisition

               Note:  The 1/(2j) refocusing period prior to acquisition is
                      always present

      phase = 1,2:  hypercomplex experiment with F1 quadrature (complex F1-FT)
                3:  TPPI experiment with F1 quadrature (complex F1-FT)
       null = nulling time for protons not attached to observed heteronucleus


   Note:  The single-bond correlation experiment can be presented
          in a phase-sensitive manner simply due to the fact that
          the total time for proton-proton scalar coupling evolution
          (2*bird + d2) is quite short. 


   s.l. patt       1986-12-31
   revised (s.f.)  1988-01-27
   revised         1988-08-01
   revised	   1989-01-24
   revised         1990-12-10 (removed assign(zero,v14))
   revised         1991-01-09 (added wfg shaped pulses) 
   revised (g.g)   1991-09-19 (made selective experiment)
   revised (r.k.)  2003-07-01 (fixed for VNMR 6.1 PSG)
*/



#include <standard.h>


pulsesequence()
{
/* VARIABLE DECLARATION */
   double at = getval("at"),
          pwxlvl = getval("pwxlvl"),
          pwx = getval("pwx"),
          j = getval("j"),
          jeff = getval("jeff"),
          bird, birdeff = 0.0,
          total_time,
          decup_time,
          null = getval("null"),
          satfrq = getval("satfrq"),
          satdly = getval("satdly"),
          satpwr = getval("satpwr"),
          selpwr = getval("selpwr"),
          selpw = getval("selpw");
   int    decup_ok;
   char   satmode[MAXSTR],
          selshape[MAXSTR],
          sspul[MAXSTR],
          satshape[MAXSTR];

   getstr("sspul", sspul);
   getstr("satmode", satmode);
   getstr("satshape", satshape);
   getstr("selshape", selshape);

   /* INITIALIZE VARIABLES */
   if (j > 0.0)
   {
      bird = 1.0 / (2.0*j);
      birdeff= 1.0 / (2.0*jeff);
   }
   else
   {
      bird = 0.0;
   }

   /* CHECK CONDITIONS */
   /* Check for correct system hardware configuration */
   if (selpw < 100e-6)
   {
      printf("selpw should be in microseconds.\n");
      abort(1);
   }
   if (!newtransamp)
   {
      printf("This sequence requires linear amplifiers on XMTR.\n");
      abort(1);
   }

   /* Check for correct DM settings */
   if ((dm[A] == 'y') || (dm[B] == 'y'))
   {
      printf("DM must be set to either 'nny' or 'nnn'.\n");
      abort(1);
   }

   /* Check for correct decoupler duty cycle */
   if (null != 0.0)
   {
      decup_time = 6.0*pwx;
   }
   else
   {
      decup_time = 2.0*pwx;
   }
   if (dm[C] == 'y')
      decup_time = decup_time + at;

   total_time = satdly + d1 + at + 2.0*bird + null;
   if (null != 0.0)
      total_time = total_time + 2.0*bird;

   decup_ok = (((decup_time / total_time) < 0.2)||(dpwr < 54.0));
   if (!decup_ok)
   {
      printf("Decoupler duty cycle must be under 20%%.\n");
      abort(1);
   }


   /* STEADY-STATE PHASECYCLING
   /* This section determines if the phase calculations trigger off of
      (SS - SSCTR) or off of CT */
   sub(ct, ssctr, v11);
   hlv(v11, v1);
   dbl(v11, v3);

   /* PHASECYCLE CALCULATION */
   hlv(v1, v4);

   dbl(v4, v8);			/* suppression of 1H 180 artifacts */
   mod2(v1,v1);			/* QIS subcycle */

   hlv(v4, v5);
   hlv(v5, v5);
   dbl(v5, v5);			/* suppression of artifacts from 2nd X 90 */

   add(v1, one, v9);
   add(two, v1, v2);		/* for composite 180 */

   add(v5, v3, oph);
   add(v8, oph, oph);		/* receiver phasecycle */

   add(oph, v1, oph);
   add(v3, v1, v3);
   add(v4, v1, v4);
   add(v5, v1, v5);		/* addition of QIS subcycle */

   if (phase1 == 2)
      incr(v3);
   if (phase1 == 3)		/* TPPI with complex FT */
      add(id2, v3, v3);

   if ((phase1 == 1) || (phase1 == 2))
   {
      dbl(id2,v14);
      add(v3,v14,v3); add(oph,v14,oph);
   }


   /* BEGIN ACTUAL PULSE SEQUENCE CODE */
   status(A);
      obspower(tpwr);
      decpower(pwxlvl);
      if (sspul[A]=='y')
      {
         rgpulse(200*pw,zero,rof1,0.0);
         rgpulse(200*pw,one,0.0,rof1);
      }
      hsdelay(d1);

      /* selective saturation period */
      if (satmode[A] == 'y')
      {
         obsoffset(satfrq);
         obspower(satpwr);
         shaped_pulse(satshape,satdly, zero, rof1, 0.0);
         obsoffset(tof);
         obspower(tpwr);
         delay(4.0e-5);
      }

   status(B);
      /* if null is 0 eliminate bird inversion pulse */
      /* bird pulse - 180 for protons not coupled to 13C */
      if ((null != 0.0) )
      {
         rcvroff();
         rgpulse(pw, v1, rof1, 0.0);
         if (pw > pwx)
         {
            delay(bird - rof1 - 1.0e-6 - 1.5*pw - pwx);
         }
         else
         {
            delay(bird - rof1 - 1.0e-6 - 2*pwx - 0.5*pw);
         }
         decrgpulse(pwx, v9, rof1, 0.0);
         simpulse(2.0*pw, 2.0*pwx, v1, v1, 1.0e-6, 0.0);
         decphase(v9);
         decrgpulse(pwx, v9, 1.0e-6, 0.0);
         if (pw > pwx)
         {
            delay(bird - rof1 - 1.0e-6 - 1.5*pw - pwx);
         }
         else
         {
            delay(bird - rof1 - 1.0e-6 - 2*pwx - 0.5*pw);
         }
         rgpulse(pw, v2, rof1, rof2);
         rcvron();


         /* nulling time for protons not coupled to 13C */
         if (satmode[B] == 'y')
         {
	    obsoffset(satfrq);
	    obspower(satpwr);
	    shaped_pulse(satshape, null, zero, rof1, 0.0);
	    obsoffset(tof);
	    obspower(tpwr);
	    delay(4.0e-5);
         }
         else
         {
	    hsdelay(null);
         }
      }

      decpower(selpwr);
      rcvroff();
      rgpulse(pw, v1, rof1, 0.0);
      delay(birdeff);
      txphase(v4);			/* presets transmitter phase */
      decshaped_pulse(selshape,selpw, v3, rof1, rof1);
      if (d2>0.0) delay(d2/2.0 -rof1-pw-(2.0*selpw/3.1414));
      decphase(v5);
      rgpulse(2.0*pw, v4, 0.0, 0.0);
      if (d2>0.0) delay(d2/2.0 -rof1-pw-(2.0*selpw/3.1414));
      decshaped_pulse(selshape,selpw, v5, rof1, rof1);
      rcvron();
      decpower(dpwr);
      delay(birdeff);
   status(C);
}
