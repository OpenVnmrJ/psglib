/* slphmqc- HMQC experiment with SLP presaturation option
         heteronuclear multiple-quantum coherence with hypercomplex and
         TPPI implementations, steady-state transients for either all or
         just the first t1 increment, pulse+receiver phasecycling during
         steady-state transients, and presaturation using the transmitter.

         [includes SLP off-resonance presaturation option]

         includes  :  FAD-Hypercomplex
                      correction for d2 > 0.0

   ref: Summers et al., j. amer. chem. soc. 108:4285-4294 (1986)


   Parameters:

     satflg = 'yn':  presaturation during relaxation period (satdly) with xmtr
              'yy':  presaturation during both relaxation period (satdly)
                     and null period (null) with xmtr
     satfrq = presaturation frequency
     satdly = saturation time during the relaxation period
     satpwr = saturation power for all periods of presaturation with xmtr
        dof = frequency for decoupling during acquisition
         hs = 'yn':  homospoil pulse (hst) during the d1 relaxation delay
              'yy':  homospoil pulse during both d1 and nulling period (null)
     slpflg = 'y' : slp presaturation
     satshape = slpsatd (h2off is positive - tof > satfrq)
                slpsatu (h2off is negative - tof < satfrq)
         pw = 90 degree xmtr pulse length for protons (the observed nucleus)
     pwxlvl = power level for decoupler pulses
        pwx = 90 degree decoupler pulse length for the heteronucleus
          j = one-bond heteronuclear coupling constant (in Hz)
      mbond = 'n':  HMQC experiment
              'y':  HMBC experiment
       jnxh = long-range coupling constant for HMBC
      taumb = delay to allow long-range heteronuclear antiphase magnetization
              to evolve for multiple-bond correlation - calculated from jnxh
              if jnxh <> 0.0 else getval'ed from taumb
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


    KRISH   Jan 13, 1992
            May 01, 1992
*/



#include <standard.h>


pulsesequence()
{
/* VARIABLE DECLARATION */
   double          ss,
                   at,
                   pwxlvl,
                   pwx,
                   j,
                   jnxh,
                   bird,
		   total_time,
		   decup_time,
                   null,
                   satfrq,
                   satdly,
                   cycle1,
                   cycle2,
                   shape_pw,
                   h2off,
		   arraydim,
		   ni,
                   satpwr,
		   taumb;
   int             iphase,
                   decup_ok,
                   t1_ok;
   char            satflg[MAXSTR],
                   slpflg[MAXSTR],
                   satshape[MAXSTR],
		   mbond[MAXSTR],
                   sspul[MAXSTR];


/* LOAD VARIABLES */
   satfrq = getval("satfrq");
   satdly = getval("satdly");
   satpwr = getval("satpwr");
   h2off = getval("h2off");
   pwxlvl = getval("pwxlvl");
   arraydim = getval("arraydim");
   ni = getval("ni");
   at = getval("at");
   pwx = getval("pwx");
   null = getval("null");
   j = getval("j");
   jnxh = getval("jnxh");
   if (jnxh > 0)
    taumb = 1/(2*jnxh);
   else
    taumb = getval("taumb");
   ss = getval("ss");
   iphase = (int) (getval("phase") + 0.5);
   getstr("satflg", satflg);
   getstr("slpflg",slpflg);
   getstr("satshape",satshape);
   getstr("sspul",sspul);
   getstr("mbond", mbond);

   if ((slpflg[0] == 'y') && (h2off != 0.0))
   {
   shape_pw = 10/h2off;
   if (shape_pw < 0.0)
    shape_pw = -shape_pw;

   cycle1 = (satdly/(shape_pw + 15.4e-6));
   cycle1 = 2.0*(double)(int)(cycle1/2.0);
   initval(cycle1,v6);

   cycle2 = (null/(shape_pw + 15.4e-6));
   cycle2 = 2.0*(double)(int)(cycle2/2.0);
   initval(cycle2,v10);
   }

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
      initval((double) (ix - 1) / (arraydim / ni), v14);
   else
   {
    if ((iphase == 1)||(iphase == 2))
     initval(2.0*(double)(((int)(d2*getval("sw1")+0.5)%2)),v14);
    else
     assign(zero, v14);
   }




/* CHECK CONDITIONS */
/* Check for correct system hardware configuration */

   if (rfwg[0] != 'y')
    {
      if ((slpflg[0] == 'y') || (slpflg[1] == 'y'))
   {
    printf("NOTE: slpflg=y requires PPM in Channel 1.  P.S.G. aborted.\n");
    abort(1);
   }
    }

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
      printf("Decoupler duty cycle must be under 20%.\n");
      abort(1);
   }

   t1_ok = ((1/(2*(getval("sw1")))) > (pw + (2*pwx/3.1416) + 2.0e-6));

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

   add(v14, v3, v3);

   if (iphase == 4)		/* TPPI with real FT */
   {
      incr(v3);
      d2 += inc2D/2;
   }

   if ((iphase==1)||(iphase==2)) 
   { 
       add(oph,v14,oph);
   }


/* BEGIN ACTUAL PULSE SEQUENCE CODE */
   status(A);
      rlpower(pwxlvl, DODEV);

      if (sspul[0] == 'y')
      {
        rgpulse(200*pw,zero,rof1,0.0);
        rgpulse(200*pw,one,0.0,rof1);
      }

      hsdelay(d1);

/* selective saturation period */
      if (satflg[0] == 'y')
      {
         rlpower(satpwr,TODEV);
         if (slpflg[0] == 'y')
         {
          if (h2off != 0.0)
          {
          rcvroff();
          delay(rof1);
          starthardloop(v6);
           shaped_pulse(satshape,shape_pw,one,0.0,1e-6);
          endhardloop();
          delay(rof2);
          rcvron();
          }
          else
          {
          if (satfrq != tof)
           offset(satfrq, TODEV);
          shaped_pulse(satshape,satdly, one, 4.0e-5, 1.0e-5);
          if (satfrq != tof)
           offset(tof, TODEV);
          }
         }
         else
         {
          if (satfrq != tof)
           offset(satfrq, TODEV);
          rgpulse(satdly, one, 4.0e-5, 1.0e-5);
          if (satfrq != tof)
           offset(tof, TODEV);
         }
         rlpower(tpwr, TODEV);
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

         if (satflg[1] == 'y')
         {
	    rlpower(satpwr, TODEV);
            if (slpflg[1] == 'y')
            {
             if (h2off!= 0.0)
             {
             rcvroff();
             delay(rof1);
             starthardloop(v10);
              shaped_pulse(satshape,shape_pw,one,0.0,1e-6);
             endhardloop();
             delay(rof2);
             rcvron();
             }
             else
             {
             if (satfrq != tof)
              offset(satfrq, TODEV);
             shaped_pulse(satshape,null, one, 4.0e-5, 1.0e-5);
             if (satfrq != tof)
              offset(tof, TODEV);
             }
            }
            else
            {
             if (satfrq != tof)
              offset(satfrq, TODEV);
             rgpulse(null, one, 4.0e-5, 1.0e-5);
             if (satfrq != tof)
              offset(tof, TODEV);
            }
	    rlpower(tpwr, TODEV);
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

      if (!t1_ok)
      {
      if (d2 <= 2*pw)			/* allows 1H 180 to occur during */
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
            delay(rof2);
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
            delay(rof2);
         }
      }
      else
      {
         decrgpulse(pwx, v3, rof1, 0.0);
         delay(d2/2.0 - pw + 0.5e-6);
         decphase(v5);
	 rgpulse(2.0*pw, v4, 0.0, 0.0);
         delay(d2/2.0 - pw + 0.5e-6);
	 decrgpulse(pwx, v5, 0.0, rof2);
      }
      }
      else
      {
         decrgpulse(pwx, v3, rof1, 1e-6);
         if (d2 > 0.0)
          delay(d2/2.0 - pw - 2.0e-6 - (2*pwx/3.1416));
         else
          delay(d2/2.0);
         decphase(v5);
         rgpulse(2.0*pw, v4, 1.0e-6, 1.0e-6);
         if (d2 > 0.0) 
          delay(d2/2.0 - pw - 2.0e-6 - (2*pwx/3.1416)); 
         else 
          delay(d2/2.0);
         decrgpulse(pwx, v5, 1.0e-6, rof2);
      }

      rcvron();
      rlpower(dpwr, DODEV);

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
