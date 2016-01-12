
/* hmqc3rf- heteronuclear multiple-quantum coherence with hypercomplex and
         TPPI implementations, steady-state transients for either all or
         just the first t1 increment, pulse+receiver phasecycling during
         steady-state transients, and presaturation using the transmitter.

   ****************************************************************************
         USES SECOND DECOUPLER FOR X NUCLEI PULSING AND DECOUPLING
   ****************************************************************************

   ref: Summers et al., j. amer. chem. soc. 108:4285-4294 (1986)


   Parameters:

     satmode = 'yn':  presaturation during relaxation period (satdly) with xmtr
              'yy':  presaturation during both relaxation period (satdly)
                     and null period (null) with xmtr
     satfrq = presaturation frequency
     satdly = saturation time during the relaxation period
     satpwr = saturation power for all periods of presaturation with xmtr
        dof = frequency for decoupling during acquisition
         hs = 'yn':  homospoil pulse (hst) during the d1 relaxation delay
              'yy':  homospoil pulse during both d1 and nulling period (null)
         pw = 90 degree xmtr pulse length for protons (the observed nucleus)
     pwx2lvl = power level for decoupler pulses
        pwx2 = 90 degree decoupler pulse length for the heteronucleus
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


   s.l. patt       31  december    1986
   revised (s.f.)  27  january     1988
   revised          1  august      1988
   revised	   24  january	   1989
   revised         10  december    1990 (removed assign(zero,v14))  */



#include <standard.h>


pulsesequence()
{
/* VARIABLE DECLARATION */
   double          ss,
                   at,
                   pwx2lvl,
                   pwx2,
                   j,
                   bird,
		   total_time,
		   decup_time,
                   null,
                   satfrq,
                   satdly,
		   arraydim,
		   ni,
                   satpwr,
		   taumb;
   int             iphase,
                   decup_ok;
   char            satmode[MAXSTR],
		   mbond[MAXSTR];


/* LOAD VARIABLES */
   satfrq = getval("satfrq");
   satdly = getval("satdly");
   satpwr = getval("satpwr");
   pwx2lvl = getval("pwx2lvl");
   arraydim = getval("arraydim");
   ni = getval("ni");
   at = getval("at");
   pwx2 = getval("pwx2");
   null = getval("null");
   j = getval("j");
   ss = getval("ss");
   iphase = (int) (getval("phase") + 0.5);
   getstr("satmode", satmode);
   getstr("mbond", mbond);

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

   if (iphase == 3)
      initval((double)((int)((ix-1)/(arraydim/ni)+1e-6)),v14); 

   initval(pwx2lvl, v10);
   initval(satpwr, v7);
   initval(tpwr, v13);
   initval(dpwr2, v6);


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
   if (dm[C] == 'y')
      decup_time = decup_time + at;

   total_time = satdly + d1 + at + 2.0*bird + null;
   if (null != 0.0)
      total_time = total_time + 2.0*bird;

   decup_ok = (((decup_time / total_time) < 0.2)||(dpwr < 45.0));
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
      power(v10, DO2DEV);
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
         dec2phase(v9);
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
      delay(bird - 0.5*pw - 0.5*pwx2 - rof1);

      if (mbond[0] == 'y')			/* one-bond J(CH)-filter */
      {
         dec2rgpulse(pwx2, v11, rof1, 0.0);
         delay(taumb - bird - rof1 - pwx2);
      }

      txphase(v4);			/* presets transmitter phase */
      if (d2 <= 2*pw)			/* allows 1H 180 to occur during */
      {					/* t1 evolution period */
         dec2phase(v3);			/* presets decoupler phase */
         delay(rof1);
         if (pwx2 >= (pw - 0.5e-6))
         {
            dec2on();
            delay(pwx2 - (2*pw - d2)/2);
            xmtron();
            delay(pw - d2/2);
            dec2off();
            dec2phase(v5);
            delay(d2);
            dec2on();
            delay(pw - d2/2);
            xmtroff();
            delay(pwx2 - (2*pw - d2)/2);
            dec2off();
         }
         else
         {
            xmtron();
            delay((2*pw - d2)/2 - pwx2);
            dec2on();
            delay(pwx2);
            dec2off();;
            dec2phase(v5);
            delay(d2);
            dec2on();
            delay(pwx2);
            dec2off();
            delay((2*pw - d2)/2 - pwx2);
            xmtroff();
         }
      }
      else
      {
         dec2rgpulse(pwx2, v3, rof1, 0.0);
         delay(d2/2.0 - pw -2*pwx2/3.1414);
         dec2phase(v5);
	 rgpulse(2.0*pw, v4, 0.0, 0.0);
         delay(d2/2.0 - pw -2*pwx2/3.1414);
	 dec2rgpulse(pwx2, v5, 0.0, 0.0);
      }

      rcvron();
      delay(rof2);
      power(v6, DO2DEV);

      if (mbond[0] == 'n')
      {
         if (pwx2 > pw)
         {
            delay(bird - 4.2e-6 - 0.5*pwx2 - rof2);
         }
         else
         {
            delay(bird - 4.2e-6 - 0.5*pw - rof2);
         }
      }

   status(C);
}
