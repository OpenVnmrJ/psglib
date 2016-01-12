/* hsqc3rf_C13dec    HSQC experiment with optional C13 decoupling in t1

              Uses decoupler for C13 pulses and decoupling
              Uses second decoupler for X pulses and decoupling
   Parameters:

      sspul = 'y':  selects for sstrim(x)-sstrim(y) sequence at the start
		    of the pulse sequence
              'n':  normal experiment
        fad = 'y':  TPPI axial-peak displacement
              'n':  standard phasecycle
      f1180 = 'y':  the first t1 point is sampled at half the t1 dwell time
              'n':  the first t1 point is sampled at t1 = 0
     satmode = 'yn':  presaturation during relaxation period (satdly) with xmtr
              'nn':  no presaturation during relaxation period (satdly)
	      'ny':  presaturation during only the null period
     satfrq = presaturation frequency
     satdly = saturation time during the relaxation period
     satpwr = saturation power for all periods of presaturation with xmtr
         hs = 'yn':  homospoil pulse (hst) during the d1 relaxation delay
       null = delay associated with the BIRD nulling
       tpwr = power level for 1H transmitter pulses
         pw = 90 degree xmtr pulse length for protons (the observed nucleus)
     pwx2lvl = power level for X decoupler pulses
        pwx2 = 90 degree decoupler pulse length for X
        jxh = one-bond heteronuclear coupling constant to X (in Hz)
    deltaxh = 1/(4*jxh) if jxh = 0.0; otherwise, the entered value is
              used; the delay used in the REVINEPT subsequences
        dm2 = 'nnn':  no broadband decoupling of X during acquisition
              'nny':  broadband heteronuclear decoupling of X during
		        acquisition
        dm  = 'nyn'   status-controlled carbonyl carbon decoupling used

     pwca   = Calpha decoupling 180 degree pulse
     pwcalvl= Calpha decoupling power level
     shape  = shapelib entry for SLP-Ca pulse shape

        dof2 = frequency for X-decoupling during acquisition
        dof  = set in the CO region of spectrum
         pw = 90 degree xmtr pulse length for protons (the observed nucleus)

      permits C13 decoupling during t1 by status control on carbonyls and
      slp-180 degree pulse on alpha carbons. A file must exist in shapelib
      to do the alpha 180. Use "convolute" to make this file, for example:

        convolute('120hard','c13dec',24,-18000), where 120 hard is a shapelib
        file that is 120 lines of 0.0  1023.0 (5 steps per microsec x 24 us =
        120 steps). 'c13dec' is the name of the shifted frequency pulse used
        by the simshaped_pulse statement, and -18000 is the effective frequency
        shift of the pulse. 

      phase = 1,2:  hypercomplex experiment with F1 quadrature (complex F1-FT)
*/


#include <standard.h>
#include <math.h>

#define MIN_J		0.1		/* Hz  */
#define MIN_DELAY	0.2e-6		/* sec */
#define MIN_NULL	0.0001		/* sec */
#define MAX_SSTRIM	0.1		/* sec */

static int	phs1[4]  = {1,1,3,3},
		phs2[2]  = {0,2},
		phs3[8]  = {0,0,0,0,2,2,2,2},
		phs4[16] = {0,0,0,0,0,0,0,0,2,2,2,2,2,2,2,2},
		phs5[16] = {0,2,2,0,2,0,0,2,2,0,0,2,0,2,2,0};

static double	d2_init = 0.0;



/*-------------------------------
|				|
|	pulsesequence()/0	|
|				|
+------------------------------*/
pulsesequence()
{
/* VARIABLE DECLARATION */
   char         satmode[MAXSTR],
		sspul[MAXSTR],
		fad[MAXSTR],shape[MAXSTR],
		f1180[MAXSTR];
   int          phase,
		satmove,
		t1_counter;
   double       ss,
		sstrim,
		sw1,
		t1evol,
                pwx2lvl,
                pwx2,pwca,pwcalvl,
                jxh,
		deltaxh,
		bird,
                satfrq,
                satdly,
                satpwr,
		null;


/* Load variables */
   satfrq = getval("satfrq");
   satdly = getval("satdly");
   satpwr = getval("satpwr");
   pwx2lvl = getval("pwx2lvl");
   pwx2 = getval("pwx2");
   jxh = getval("jxh");
   deltaxh = getval("deltaxh");
   ss = getval("ss");
   sw1 = getval("sw1");
   sstrim = getval("sstrim");
   null = getval("null");
   phase = (int) (getval("phase") + 0.5);

   getstr("sspul", sspul);
   getstr("f1180", f1180);
   getstr("fad", fad);
   getstr("satmode", satmode);

    pwca = getval("pwca");
  pwcalvl= getval("pwcalvl");
  getstr("shape",shape);

/* Load phase tables */
   settable(t1, 4, phs1);
   settable(t2, 2, phs2);
   settable(t3, 8, phs3);
   settable(t4, 16, phs4);
   settable(t5, 16, phs5);



/* Adjust delays */
   if (jxh > MIN_J)
   {
      bird = 1/(2*jxh);
      deltaxh = 0.4*bird;
   }
   else
   {
      bird = 0.0;
   }


/* Check for 1H frequency change */
   satmove = ( fabs(tof - satfrq) >= 0.1 );


/* Check for correct 'dm/dm2' and 'dpwr/dpwr2' settings */
      if ( (dm[A] == 'y') || (dm[C] == 'y') )
      {
         text_error("DM must be set to either 'nnn' or 'nyn'\n");
         abort(1);
      }
      if ( (dm2[A] == 'y') || (dm2[B] == 'y'))
      { 
         text_error("DM2 must be set to either 'nnn' or 'nny'\n");
         abort(1);
      }

   if((dpwr2 > 50) || (dpwr>40))
   {
      text_error("dpwr2 or dpwr is too large, reset\n");
      abort(1);
   }

   if(sstrim > MAX_SSTRIM)
   {
      text_error("`sstrim` is > maximum value\n");
      abort(1);
   }


/* Determine steady-state mode */
   if (ss < 0)
   {
      ss *= (-1);
      initval(ss, ssval);
      initval(ss, ssctr);
   }


/* Phase incrementation for hypercomplex 2D data */
   if (phase == 2)
      tsadd(t2, 1, 4);


/* FAD phase incrementation */
   if (fad[A] == 'y')
   {
      if (ix == 1)
        d2_init = d2;
      t1_counter = (int) ( (d2 - d2_init)*sw1 + 0.5 );
      if (t1_counter % 2)
      {
         tsadd(t2, 2, 4);
         tsadd(t5, 2, 4);	/* receiver phase cycle */
      }
   }


/* BEGIN ACTUAL PULSE SEQUENCE CODE */
   status(A);
      rlpower(tpwr, TODEV);
      rlpower(pwx2lvl, DO2DEV);

      if (sspul[A] == 'y')
      {
         rlpower(tpwr-12, TODEV);
         rgpulse(sstrim, zero, rof1, 1.0e-6);
         rgpulse(sstrim, one, rof1, rof2);
         rlpower(tpwr, TODEV);
      }
         hsdelay(d1);

/* selective saturation period */
      if (satmode[A] == 'y')
      {
         if (satmove)
            offset(satfrq, TODEV);

         rlpower(satpwr, TODEV);
         rgpulse(satdly, zero, 4.0e-5, 0.2e-6);
         if (satmove)
            offset(tof, TODEV);

         rlpower(tpwr, TODEV);
         delay(1.0e-5);
      }


/* Bird pulse and nulling period for both C13 and N15 */
      if (null > MIN_NULL)
      {
         rgpulse(pw, zero, rof1, 0.0);
         delay(bird - rof1 - 1.0e-6 - 2*pwx2 - 0.5*pw);

         rcvroff();
         dec2rgpulse(pwx2, one,0.0,0.0);
         sim3pulse(2*pw,0.0, 2*pwx2, zero,zero, zero, 1.0e-6, 0.0);
         dec2rgpulse(pwx2, one,0.0,0.0);
         rcvron();

         delay(bird - rof1 - 1.0e-6 - 2*pwx2 - 0.5*pw);
         rgpulse(pw, two, rof1, 1.0e-6);
         
         
         if (satmode[B] == 'y')
         {
            if (satmove)
               offset(satfrq, TODEV);

            rlpower(satpwr, TODEV);
            rgpulse(null, zero, 1.0e-5, 0.2e-6);
            if (satmove)
               offset(tof, TODEV);

            rlpower(tpwr, TODEV);
            delay(1.0e-5);
         }
         else
         {
            delay(null);
         }
      }


      rcvroff();
      rgpulse(pw, zero, rof1, 0.0);
      delay(deltaxh - pwx2);
      sim3pulse(2*pw,0.0, 2*pwx2, zero,zero, zero, 0.0, 0.0);
      delay(deltaxh - pwx2);
      sim3pulse(pw, 0.0,pwx2, t1,zero, t2, 0.0, 0.0);

/* Calculate t1 delay */
      t1evol = d2;
      if (f1180[A] == 'y')
         t1evol += 0.5/sw1;


      if (t1evol > MIN_DELAY)
      {
         t1evol -= pwca + (4*pwx2/M_PI);
         if (t1evol < MIN_DELAY)
            t1evol = 0.0;
      }
    if (d2>(pwca +53.2e-6 +8.4e-6))
     {
     status(B);                                 /* turn on CO decoupling */
         delay( t1evol/2 -26.6e-6 );
     status(A);                                 /* turn off CO decoupling */
         rlpower(pwcalvl,DODEV);
         simshaped_pulse("hard",shape,2.0*pw,pwca,zero,zero,0.0,0.0);
         rlpower(dpwr,DODEV);
     status(B);                                 /* turn on CO decoupling */
         delay(t1evol/2 -26.6e-6);
     status(A);
     }
      sim3pulse(pw, 0.0,pwx2, t3, zero,t4, 0.0, 0.0);
      delay( deltaxh - pwx2 - (2*pwx2/M_PI) );
      sim3pulse(2*pw,0.0, 2*pwx2, zero,zero, zero, 0.0, 0.0);

      rlpower(dpwr2, DO2DEV);		/* X decoupling power level */
      delay(rof2);
      rcvron();
      delay(deltaxh - pwx2 - POWER_DELAY - rof2);
      setreceiver(t5);
     status(C);
}
