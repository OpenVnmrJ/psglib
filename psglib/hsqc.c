/* overbdn1.c  -  heteronuclear Overbodenhausen experiment using REVINEPT


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
     pwxlvl = power level for X decoupler pulses
        pwx = 90 degree decoupler pulse length for X
     xrfdev = RF device for the X heteronucleus(2 for 1st dec., 3 for 2nd dec.)
        jxh = one-bond heteronuclear coupling constant to X (in Hz)
    deltaxh = 1/(4*jxh) if jxh = 0.0; otherwise, the entered value is
              used; the delay used in the REVINEPT subsequences
    dm(dm2) = 'nnnnn':  no broadband decoupling of X during acquisition
              'nnnny':  broadband heteronuclear decoupling of X during
		        acquisition


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
		fad[MAXSTR],
		f1180[MAXSTR];
   int          phase,
		satmove,
		t1_counter,
		xrfdev;
   double       ss,
		sstrim,
		sw1,
		t1evol,
		xdecpwr = 0.0,	/* safety precaution */
                pwxlvl,
                pwx,
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
   pwxlvl = getval("pwxlvl");
   pwx = getval("pwx");
   jxh = getval("jxh");
   deltaxh = getval("deltaxh");
   ss = getval("ss");
   sw1 = getval("sw1");
   sstrim = getval("sstrim");
   null = getval("null");
   phase = (int) (getval("phase") + 0.5);
   xrfdev = (int) (getval("xrfdev") + 0.5);

   getstr("sspul", sspul);
   getstr("f1180", f1180);
   getstr("fad", fad);
   getstr("satmode", satmode);


/* Load phase tables */
   settable(t1, 4, phs1);
   settable(t2, 2, phs2);
   settable(t3, 8, phs3);
   settable(t4, 16, phs4);
   settable(t5, 16, phs5);


/* Check X RF device */
   if ( (xrfdev != DODEV) && (xrfdev != DO2DEV) )
   {
      text_error("invalid X RF device\n");
      abort(1);
   }


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
   if (xrfdev == DODEV)
   {
      if ( (dm[A] == 'y') || (dm[B] == 'y') || (dm[C] == 'y') ||
		(dm[D] == 'y') )
      {
         text_error("DM must be set to either 'nnnnn' or 'nnnny'\n");
         abort(1);
      }

      xdecpwr = dpwr;
   }
   else
   {
      if ( (dm2[A] == 'y') || (dm2[B] == 'y') || (dm2[C] == 'y') ||
		(dm2[D] == 'y') )
      { 
         text_error("DM2 must be set to either 'nnnnn' or 'nnnny'\n");
         abort(1);
      }

      xdecpwr = dpwr2;
   }


   if (sstrim > MAX_SSTRIM)
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
      rlpower(pwxlvl, xrfdev);

      if (sspul[A] == 'y')
      {
         rgpulse(sstrim, zero, rof1, 1.0e-6);
         rgpulse(sstrim, one, rof1, rof2);
         hsdelay(d1);
      }
      else
      {
         hsdelay(d1);
      }

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


   status(B);
/* Bird pulse and nulling period for both C13 and N15 */
      if (null > MIN_NULL)
      {
         rgpulse(pw, zero, rof1, 0.0);
         delay(bird - rof1 - 1.0e-6 - 2*pwx - 0.5*pw);

         rcvroff();
         genrgpulse(pwx, one, rof1, 0.0, xrfdev);
         gensim2pulse(2*pw, 2*pwx, zero, zero, 1.0e-6, 0.0, TODEV, xrfdev);
         genrgpulse(pwx, one, 1.0e-6, 0.0, xrfdev);
         rcvron();

         delay(bird - rof1 - 1.0e-6 - 2*pwx - 0.5*pw);
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


   status(C);
      rcvroff();
      rgpulse(pw, zero, rof1, 0.0);

      genqdphase(zero, xrfdev);
      txphase(zero);
      delay(deltaxh - pwx);
      gensim2pulse(2*pw, 2*pwx, zero, zero, 0.0, 0.0, TODEV, xrfdev);
      txphase(t1);
      genqdphase(t2, xrfdev);
      delay(deltaxh - pwx);
      gensim2pulse(pw, pwx, t1, t2, 0.0, 0.0, TODEV, xrfdev);

      txphase(zero);
      genqdphase(t4, xrfdev);

/* Calculate t1 delay */
      t1evol = d2;
      if (f1180[A] == 'y')
         t1evol += 0.5/sw1;


      if (t1evol > MIN_DELAY)
      {
         t1evol -= 2*pw + (4*pwx/M_PI);
         if (t1evol < MIN_DELAY)
            t1evol = 0.0;
      }

      delay(t1evol/2);
      rgpulse(2*pw, zero, 0.0, 0.0);
      txphase(t3);
      delay(t1evol/2);


   status(D);
      gensim2pulse(pw, pwx, t3, t4, 0.0, 0.0, TODEV, xrfdev);
      txphase(zero);
      genqdphase(zero, xrfdev);
      delay( deltaxh - pwx - (2*pwx/M_PI) );
      gensim2pulse(2*pw, 2*pwx, zero, zero, 0.0, 0.0, TODEV, xrfdev);

      rlpower(xdecpwr, xrfdev);		/* X decoupling power level */
      delay(rof2);
      rcvron();
      delay(deltaxh - pwx - POWER_DELAY - rof2);


   status(E);
      setreceiver(t5);
}
