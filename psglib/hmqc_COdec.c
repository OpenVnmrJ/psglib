
/* hmqc_COdec- heteronuclear multiple-quantum coherence with hypercomplex and
         TPPI implementations, steady-state transients for either all or
         just the first t1 increment, pulse+receiver phasecycling during
         steady-state transients, and presaturation using the transmitter.

         uses status-controlled n15 decoupling during t1 and an SLP-180
         degree band-selective pulse on the carbonyls in t1 for decoupling.
         Must have a shapelib entry for "shape" which is frequency-shifted 
         enough from dof to be centered in the carbonyl carbon region.

   ref: Summers et al., j. amer. chem. soc. 108:4285-4294 (1986)


   Parameters:
      sspul = 'y' does trim(x)-trim(y) before d1(recommended for hmbc)
     satmode = 'yn':  presaturation during relaxation period (satdly) with xmtr
              'yy':  presaturation during both relaxation period (satdly)
                     and null period (null) with xmtr
     satfrq = presaturation frequency
     satdly = saturation time during the relaxation period
     satpwr = saturation power for all periods of presaturation with xmtr
     shape  = shapelib file for SLP-carbonyl decoupling pulse
    pwcolvl = power level for carbonyl 180
      pwco  = pulse length for carbonyl 180
        dof = frequency for decoupling during acquisition
         hs = 'yn':  homospoil pulse (hst) during the d1 relaxation delay
              'yy':  homospoil pulse during both d1 and nulling period (null)
         pw = 90 degree xmtr pulse length for protons (the observed nucleus)
     pwxlvl = power level for decoupler pulses
        pwx = 90 degree decoupler pulse length for the heteronucleus
          j = one-bond heteronuclear coupling constant (in Hz)
      mbond = 'n':  HMQC experiment
              'y':  HMBC experiment
      taumb = delay to allow long-range heteronuclear antiphase magnetization
              to evolve for multiple-bond correlation
         dm = 'nnn':  no broadband decoupling of heteronucleus during
                      acquisition; this is mandatory if mbond='y'
              'nny':  broadband heteronuclear decoupling during acquisition
                      (recommended for mbond='n')

               Note:  The 1/(2j) refocusing period prior to acquisition is
                      always present, regardless of dm, unless mbond='y'

      phase = 1,2:  hypercomplex experiment with F1 quadrature (complex F1-FT)
                3:  TPPI experiment with F1 quadrature (complex F1-FT)
              1,4:  TPPI experiment with F1 quadrature (real F1-FT)
       null = nulling time for protons not attached to observed heteronucleus
         nt = multiple of  4 (minimum)
	      multiple of 32 (maximum for mbond='n')
              multiple of 64 (maximum for mbond='y')
	      multiple of 16 (minimum recommended for mbond='n')
              multiple of 32 (minimum recommended for mbond='y')


   Note:  The single-bond correlation experiment can be presented
          in a phase-sensitive manner simply due to the fact that
          the total time for proton-proton scalar coupling evolution
          (2*bird + d2) is quite short.  The multiple-bond correla-
          tion experiment requires an absolute value display in F2
          but not in F1.


    added sspul and d2 correction to vnmr 4.1 version of hmqc.c  */



#include <standard.h>


pulsesequence()
{
/* VARIABLE DECLARATION */
   double          ss,
                   at,
                   pwxlvl,
                   pwx,
                   j,
                   bird,
		   total_time,
		   decup_time,
                   null,
                   satfrq,
                   satdly,
		   arraydim,
		   ni,
                   pwcolvl,
                   pwco,
                   satpwr,
		   taumb;
   int             iphase,
                   decup_ok;
   char            satmode[MAXSTR],
                   shape[MAXSTR],
                   sspul[MAXSTR],
		   mbond[MAXSTR];


/* LOAD VARIABLES */
   satfrq = getval("satfrq");
   satdly = getval("satdly");
   satpwr = getval("satpwr");
   pwxlvl = getval("pwxlvl");
   arraydim = getval("arraydim");
   pwco   = getval("pwco");
   pwcolvl=getval("pwcolvl");
   ni = getval("ni");
   at = getval("at");
   pwx = getval("pwx");
   null = getval("null");
   j = getval("j");
   ss = getval("ss");
   iphase = (int) (getval("phase") + 0.5);
   getstr("satmode", satmode);
   getstr("sspul", sspul);
   getstr("mbond", mbond);
   getstr("shape", shape);

   if (mbond[0] == 'y')
      taumb = getval("taumb");


/* INITIALIZE VARIABLES */
   if (j > 0.0)
   {
      bird = 1.0 / (2.0*j);
   }
   else
   {
      bird = 0.0;
   }

/* For TPPI,  initialize v14 with the t1 increment number */
   if (iphase == 3)
      initval((double) ((int) ((ix - 1) / (arraydim / ni) + 1e-6)), v14);

   initval(pwxlvl, v10);
   initval(satpwr, v7);
   initval(tpwr, v13);

   if (newdecamp)
      initval(dpwr, v6);


/* CHECK CONDITIONS */
/* Check for correct system hardware configuration */
   if (!newdec)
   {
      printf("This sequence requires direct synthesis RF on DEC.\n");
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
   if ((mbond[0] == 'y') && (dm[C] == 'y'))
   {
      printf("DM must be set to 'nnn' for multiple-bond correlation.\n");
      abort(1);
   }

/* Check for ROF1 minimum value */
   if ((rof1 < 9.9e-6) && (ix == 1))
      printf("Warning:  ROF1 is less than 10 us\n");

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

   decup_ok = (((decup_time / total_time) < 0.2) && (dpwr < 55.0));
   if (!decup_ok)
   {
      printf("Decoupler duty cycle must be under 20%%.\n");
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
      sub(ssval, ssctr, v12);	/* v12 = 0,...,ss-1 */
      hlv(v12, v1);
      dbl(v12, v3);
   endif(ssctr);


/* PHASECYCLE CALCULATION */
   if (mbond[0] == 'y')
   {
      hlv(v1, v11);
      hlv(v11, v4);
      dbl(v11, v11);		/* J-filter phase subcycle */
   }
   else
   {
      hlv(v1, v4);
   }

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


/* BEGIN ACTUAL PULSE SEQUENCE CODE */
   status(A);
      if (newdecamp)
      {
         power(v10, DODEV);
      }
      else
      {
         declvlon();
      }
      if (sspul[A] == 'y')
      {
       power(v13,TODEV);
       rgpulse(200*pw,zero,rof1,0.0);
       rgpulse(200*pw,one,0.0,rof1);
       }
      hsdelay(d1);

/* selective saturation period */
      if (satmode[A] == 'y')
      {
         offset(satfrq, TODEV);
         power(v7, TODEV);
         rgpulse(satdly, zero, 4.0e-5, 1.0e-5);
         offset(tof, TODEV);
         power(v13, TODEV);
         delay(4.0e-5);
      }


/* if null is 0 eliminate bird inversion pulse */
/* bird pulse - 180 for protons not coupled to 13C */

   status(B);
      if ((null != 0.0) && (mbond[0] == 'n'))
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
	    offset(satfrq, TODEV);
	    power(v7, TODEV);
	    rgpulse(null, zero, 4.0e-5, 1.0e-5);
	    offset(tof, TODEV);
	    power(v13, TODEV);
	    delay(4.0e-5);
         }
         else
         {
	    hsdelay(null);
         }
      }


      rcvroff();
      rgpulse(pw, v1, rof1, 0.0);
      delay(bird - 0.5*pw - 0.5*pwx - rof1);

      if (mbond[0] == 'y')			/* one-bond J(CH)-filter */
      {
         decrgpulse(pwx, v11, rof1, 0.0);
         delay(taumb - bird - rof1 - pwx);
      }

      txphase(v4);			/* presets transmitter phase */
      if (d2 <= pwco)			/* allows 1H 180 to occur during */
      {					/* t1 evolution period */
         decphase(v3);			/* presets decoupler phase */
         delay(rof1);
         if (pwx >= (pw - 0.5e-6))
         {
            decon();
            delay(pwx - (2*pw - d2)/2);
            xmtron();
            delay(pw - d2/2);
            decoff();
            decphase(v5);
            delay(d2);
            decon();
            delay(pw - d2/2);
            xmtroff();
            delay(pwx - (2*pw - d2)/2);
            decoff();
         }
         else
         {
            xmtron();
            delay((2*pw - d2)/2 - pwx);
            decon();
            delay(pwx);
            decoff();
            decphase(v5);
            delay(d2);
            decon();
            delay(pwx);
            decoff();
            delay((2*pw - d2)/2 - pwx);
            xmtroff();
         }
      }
      else
      {
           decrgpulse(pwx, v3, rof1, 0.0);
           delay(d2/2.0  -pwco/2 -18.4e-6+ 0.5e-6 -2.0*pwx/3.1416);
           rlpower(pwcolvl,DODEV);
           decphase(v5);
            simshaped_pulse("hard",shape,2.0*pw,pwco,v4,zero,0.0,0.0);
           rlpower(pwxlvl,DODEV);
           delay(d2/2.0  -pwco/2 -18.4e-6+ 0.5e-6 -2.0*pwx/3.1416);
           decrgpulse(pwx, v5, 0.0, 0.0);
      }

      rcvron();
      delay(rof2);
      if (newdecamp)
      {
         power(v6, DODEV);
      }
      else
      {
         declvloff();
      }

      if (mbond[0] == 'n')
      {
         if (pwx > pw)
         {
            delay(bird - 4.2e-6 - 0.5*pwx - rof2);
         }
         else
         {
            delay(bird - 4.2e-6 - 0.5*pw - rof2);
         }
      }

   status(C);
}
