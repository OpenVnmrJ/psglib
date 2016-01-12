/* selhmqc3rf- selective HMQC with hypercomplex and
         TPPI implementations, steady-state transients for either all or
         just the first t1 increment, pulse+receiver phasecycling during
         steady-state transients, and presaturation using the transmitter
         using the waveform generator

           USES THE THIRD CHANNEL FOR X PULSES AND DECOUPLING

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
     selpw= selective decoupler pulse length (in usec) at dof2
   selshape = selective decoupler (wfg) shape
        do2f = frequency for decoupling during acquisition
         hs = 'yn':  homospoil pulse (hst) during the d1 relaxation delay
              'yy':  homospoil pulse during both d1 and nulling period (null)
         pw = 90 degree xmtr pulse length for protons (the observed nucleus)
     pwx2lvl = power level for decoupler pulses
        pwx2 = 90 degree decoupler pulse length for the heteronucleus
          j = one-bond heteronuclear coupling constant (in Hz)
       jeff = effective j (longer than j to account for precession during selpw)
         dm2 = 'nnn':  no broadband decoupling of heteronucleus during
                      acquisition; 
              'nny':  broadband heteronuclear decoupling during acquisition

               Note:  The 1/(2j) refocusing period prior to acquisition is
                      always present

      phase = 1,2:  hypercomplex experiment with F1 quadrature (complex F1-FT)
                3:  TPPI experiment with F1 quadrature (complex F1-FT)
              1,4:  TPPI experiment with F1 quadrature (real F1-FT)
       null = nulling time for protons not attached to observed heteronucleus


   Note:  The single-bond correlation experiment can be presented
          in a phase-sensitive manner simply due to the fact that
          the total time for proton-proton scalar coupling evolution
          (2*bird + d2) is quite short. 


   s.l. patt       31  december    1986
   revised (s.f.)  27  january     1988
   revised          1  august      1988
   revised	   24  january	   1989
   revised         10  december    1990 (removed assign(zero,v14))
   revised          9  january     1991 (added wfg shaped pulses) 
   revised (g.g)   19  september   1991 (made selective experiment) */



#include <standard.h>
#include <math.h>


pulsesequence()
{
/* VARIABLE DECLARATION */
   double          ss,
                   at,
                   pwx2lvl,
                   pwx2,
                   j , jeff,
                   bird,birdeff,
		   total_time,
		   decup_time,
                   null,
                   satfrq,
                   satdly,
		   arraydim,
		   ni,
                   satpwr,
                   selpwr,
                   selpw;
   int             iphase,
                   decup_ok;
   char            satmode[MAXSTR],
                   selshape[MAXSTR],
                   sspul[MAXSTR],
                   satshape[MAXSTR];


/* LOAD VARIABLES */
   satfrq = getval("satfrq");
   satdly = getval("satdly");
   satpwr = getval("satpwr");
   selpwr = getval("selpwr");
   selpw= getval("selpw");
   pwx2lvl = getval("pwx2lvl");
   arraydim = getval("arraydim");
   ni = getval("ni");
   at = getval("at");
   pwx2 = getval("pwx2");
   null = getval("null");
   j = getval("j");
   jeff = getval("jeff");
   ss = getval("ss");
   iphase = (int) (getval("phase") + 0.5);
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

   if (iphase == 3)
      initval((double)((int)((ix-1)/(arraydim/ni)+1.0e-6)), v14);

   if (newdecamp)
      initval(dpwr2, v6);


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

/* Check for correct DM2 and DM settings */
   if ((dm2[A] == 'y') || (dm2[B] == 'y') || (dm[A] =='y') || (dm[B] == 'y'))
   {
      printf("DM2 and DM must be set to either 'nny' or 'nnn'.\n");
      abort(1);
   }

/* Check for ROF1 minimum value */
   if ((rof1 < 9.9e-6) && (ix == 1))
      printf("Warning:  ROF1 is less than 10 us\n");

/* Check for correct decoupler duty cycle */
   if (null != 0.0)
   {
      decup_time = 6.0*pwx2;
   }
   else
   {
      decup_time = 2.0*pwx2;
   }
   if (dm2[C] == 'y')
      decup_time = decup_time + at;

   total_time = satdly + d1 + at + 2.0*bird + null;
   if (null != 0.0)
      total_time = total_time + 2.0*bird;

   decup_ok = (((decup_time / total_time) < 0.2)||(dpwr2 < 54.0));
   if (!decup_ok)
   {
      printf("Decoupler duty cycle must be under 20%.\n");
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
   initval(ss, ssval);
   initval(ss, ssctr);


/* STEADY-STATE PHASECYCLING
/* This section determines if the phase calculations trigger off of (SS - SSCTR)
   or off of CT */

   ifzero(ssctr);
      hlv(ct, v1);
      dbl(ct, v3);		/* suppression of artifacts from 1st X 90 */
   elsenz(ssctr);
      sub(ssval, ssctr, v11);	/* v11 = 0,...,ss-1 */
      hlv(v11, v1);
      dbl(v11, v3);
   endif(ssctr);


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

   if (iphase == 2)
      incr(v3);
   if (iphase == 3)		/* TPPI with complex FT */
      add(v14, v3, v3);

   if ((iphase == 1) || (iphase == 2))
   {
      initval(2.0*(double)((int)(d2*getval("sw1")+0.5)%2),v14);
      add(v3,v14,v3); add(oph,v14,oph);
   }
   if (iphase == 4)		/* TPPI with real FT */
   {
      incr(v3);
      d2 += inc2D/2;
   }


   initval(pwx2lvl, v10);
   initval(satpwr, v7);
   initval(tpwr, v13);
   initval(selpwr,v12);
/* BEGIN ACTUAL PULSE SEQUENCE CODE */
   status(A);
      power(v13,TODEV);
      power(v10, DO2DEV);
     if (sspul[A]=='y')
      {
        rgpulse(200*pw,zero,rof1,0.0);
        rgpulse(200*pw,one,0.0,rof1);
      }
      hsdelay(d1);

/* selective saturation period */
      if (satmode[A] == 'y')
      {
         offset(satfrq, TODEV);
         power(v7, TODEV);
         shaped_pulse(satshape,satdly, zero, 4.0e-5, 1.0e-5);
         offset(tof, TODEV);
         power(v13, TODEV);
         delay(4.0e-5);
      }


/* if null is 0 eliminate bird inversion pulse */
/* bird pulse - 180 for protons not coupled to 13C */

   status(B);
      if ((null != 0.0) )
      {
         rcvroff();
         rgpulse(pw, v1, rof1, 0.0);
         if (pw > pwx2)
         {
            delay(bird - rof1 - 1.0e-6 - 1.5*pw - pwx2);
         }
         else
         {
            delay(bird - rof1 - 1.0e-6 - 2*pwx2 - 0.5*pw);
         }
         dec2rgpulse(pwx2, v9, rof1, 0.0);
         sim3pulse(2.0*pw, 0.0, 2.0*pwx2, v1, zero, v1, 1.0e-6, 0.0);
         decphase(v9);
         dec2rgpulse(pwx2, v9, 1.0e-6, 0.0);
         if (pw > pwx2)
         {
            delay(bird - rof1 - 1.0e-6 - 1.5*pw - pwx2);
         }
         else
         {
            delay(bird - rof1 - 1.0e-6 - 2*pwx2 - 0.5*pw);
         }
         rgpulse(pw, v2, rof1, rof2);
         rcvron();


/* nulling time for protons not coupled to 13C */

         if (satmode[B] == 'y')
         {
	    offset(satfrq, TODEV);
	    power(v7, TODEV);
	    shaped_pulse(satshape,null, zero, 4.0e-5, 1.0e-5);
	    offset(tof, TODEV);
	    power(v13, TODEV);
	    delay(4.0e-5);
         }
         else
         {
	    hsdelay(null);
         }
      }

      power(v12,DO2DEV);
      rcvroff();
      rgpulse(pw, v1, rof1, 0.0);
      delay(birdeff);
      txphase(v4);			/* presets transmitter phase */
      dec2shaped_pulse(selshape,selpw, v3, rof1, rof1);
      if (d2>0.0) delay(d2/2.0 -rof1-pw-(2.0*selpw/3.1414));
      decphase(v5);
      rgpulse(2.0*pw, v4, 0.0, 0.0);
      if (d2>0.0) delay(d2/2.0 -rof1-pw-(2.0*selpw/3.1414));
      dec2shaped_pulse(selshape,selpw, v5, rof1, rof1);
      rcvron();
      power(v6, DO2DEV);
      delay(birdeff);
   status(C);
}
