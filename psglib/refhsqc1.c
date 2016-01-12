/* refhsq.c      -  heteronuclear Overbodenhausen experiment using
                    partially modified refocused INEPT and broadband
		    decoupling of H1 during t1; phasecycled the H1
		    purge pulse; this sequence has RF duality.


   Parameters:

      sspul = 'y':  selects for Trim(x)-Trim(y) sequence at the start
		    of the pulse sequence
              'n':  normal experiment
        fad = 'y':  TPPI axial-peak displacement
              'n':  standard phasecycle
      f1180 = 'y':  the first t1 point is sampled at half the t1 dwell time
              'n':  the first t1 point is sampled at t1 = 0
     satflg = 'y':  presaturation during relaxation period (satdly) with xmtr
              'n':  no presaturation during relaxation period (satdly)
     satfrq = presaturation frequency
     satdly = saturation time during the relaxation period
     satpwr = saturation power for all periods of presaturation with xmtr
         hs = 'yn':  homospoil pulse (hst) during the d1 relaxation delay
       tpwr = power level for H1 transmitter pulses
         pw = 90 degree xmtr pulse length for protons (the observed nucleus)
    hdshape = H1 programmable decoupling shape
      hdpwr = power level for H1 programmable decoupling
       hd90 = 90 degree pulse length for H1 programmable decoupling shape
              at power level 'hdpwr'
      hdres = tip-angle resolution for the H1 programmable decoupling
              shape
     pwxlvl = power level for X decoupler pulses
        pwx = 90 degree decoupler pulse length for X
     xrfdev = RF device for the X heteronucleus
        jxh = one-bond heteronuclear coupling constant to X (in Hz)
    deltaxh = 1/(4*jxh) if deltaxh = 0.0; otherwise, the entered value is
              used; one of the delays used in the refocussed INEPT
              subsequences
      tauxh = 0.8/(4*jxh) if tauxh = 0.0; otherwise, the entered value is
              used; the other delay used in the refocussed INEPT
              subsequences
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


static int	phs1[2]   = {0,2},
		phs2[4]   = {1,1,3,3},
		phs3[16]  = {0,0,0,0,0,0,0,0,2,2,2,2,2,2,2,2},
		phs4[32]  = {0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
			     2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2},
		phs5[128] = {0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
			     0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
			     1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,
			     1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,
			     2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,
			     2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,
			     3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,
			     3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3},
		phs6[64]  = {0,2,2,0,0,2,2,0,2,0,0,2,2,0,0,2,
			     2,0,0,2,2,0,0,2,0,2,2,0,0,2,2,0,
			     2,0,0,2,2,0,0,2,0,2,2,0,0,2,2,0,
			     0,2,2,0,0,2,2,0,2,0,0,2,2,0,0,2},
		phs7[8]   = {1,1,1,1,3,3,3,3};

static double	d2_init = 0.0;


/*-------------------------------
|				|
|	pulsesequence()/0	|
|				|
+------------------------------*/
pulsesequence()
{
/* VARIABLE DECLARATION */
   char         satflg[MAXSTR],
		sspul[MAXSTR],
		fad[MAXSTR],
		f1180[MAXSTR],
		hdshape[MAXSTR],
		*dmx;
   int          phase,
		satmove,
		t1_counter,
		xrfdev;
   double       ss,
		sstrim,
                ni,
		sw1,
		t1evol,
                pwxlvl,
                pwx,
                jxh,
		tauxh,
		deltaxh,
		bird,
                satfrq,
                satdly,
                satpwr,
		hdpwr,
		hd90,
		hdres;


/* Load variables */
   satfrq = getval("satfrq");
   satdly = getval("satdly");
   satpwr = getval("satpwr");
   pwxlvl = getval("pwxlvl");
   pwx = getval("pwx");
   jxh = getval("jxh");
   tauxh = getval("tauxh");
   deltaxh = getval("deltaxh");
   ss = getval("ss");
   sw1 = getval("sw1");
    ni = getval("ni");
   sstrim = getval("sstrim");
   hd90 = getval("hd90");
   hdres = getval("hdres");
   hdpwr = getval("hdpwr");
   phase = (int) (getval("phase") + 0.5);
   xrfdev = (int) (getval("xrfdev") + 0.5);

   getstr("sspul", sspul);
   getstr("f1180", f1180);
   getstr("fad", fad);
   getstr("satflg", satflg);
   getstr("hdshape", hdshape);	/* H1 t1 decoupling pattern	*/


/* Load phase tables */
   settable(t1, 2, phs1);
   settable(t2, 4, phs2);
   settable(t3, 16, phs3);
   settable(t4, 32, phs4);
   settable(t5, 128, phs5);
   settable(t6, 64, phs6);
   settable(t7, 8, phs7);


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
      tauxh = bird/2;
      deltaxh = 0.4*bird;
   }
   else
   {
      bird = 0.0;
   }


/* Check for H1 frequency and C13 power changes */
   satmove = ( fabs(tof - satfrq) >= 0.1 );


/* Check for correct 'dm/dm2' and 'dpwr/dpwr2' settings */
   dmx = ( (xrfdev == DODEV) ? dm : dm2 );

   if ( (dmx[A] == 'y') || (dmx[B] == 'y') || (dmx[C] == 'y') ||
	(dmx[D] == 'y') )
   {
      text_error("dm(X) must be set to either 'nnnnn' or 'nnnny'\n");
      abort(1);
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
      tssub(t1, 1, 4);


/* FAD phase incrementation */
   if (fad[A] == 'y')
   {
      if (ix == 1)
        d2_init = d2;
      t1_counter = (int) ( (d2 - d2_init)*sw1 + 0.5 );
      if (t1_counter % 2)
      {
         tsadd(t1, 2, 4);
         tsadd(t6, 2, 4);	/* receiver phase cycle */
      }
   }


/* BEGIN ACTUAL PULSE SEQUENCE CODE */
   status(A);
      rlpower(pwxlvl, xrfdev);
      rlpower(hdpwr, TODEV);
      obsprgon(hdshape, hd90, hdres);
      xmtron();
      delay((ni/sw1)-d2);  /* total 1H decoupling time independent of d2 */
      xmtroff();
      obsprgoff();
      rlpower(tpwr, TODEV);
 
      if (sspul[A] == 'y') /* establish same inital magnetization for every d2 */
      {
         rgpulse(sstrim, zero, rof1, 1.0e-6);
         rgpulse(sstrim, one, rof1, rof2);
      }
      hsdelay(d1);
/* selective saturation period */
      if (satflg[A] == 'y')
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
/* does nothing now */


   status(C);
      rcvroff();
      rgpulse(pw, zero, rof1, 0.0);

      genqdphase(zero, xrfdev);
      delay(deltaxh - pwx);
      gensim2pulse(2*pw, 2*pwx, zero, zero, 0.0, 5.0e-6, TODEV, xrfdev);

      txphase(t2);
      genqdphase(t1, xrfdev);
      delay( deltaxh - pwx - 5.0e-6 - (pwx - pw)/2 );
      gensim2pulse(pw, pwx, t2, t1, 0.0, 5.0e-6, TODEV, xrfdev);

      genqdphase(t5, xrfdev);
      rlpower(hdpwr, TODEV);
      txphase(zero);
      delay(2*tauxh - POWER_DELAY - PRG_START_DELAY - 5.0e-6);


/* Calculate t1 delay */
      t1evol = d2;
      if (f1180[A] == 'y')
         t1evol += 0.5/sw1;

/* Start H1 decoupling */
      obsprgon(hdshape, hd90, hdres);
      xmtron();
      genrgpulse(2*pwx, t5, 0.0, 0.0, xrfdev);

      delay(t1evol);

      xmtroff();
      obsprgoff();
      rlpower(tpwr, TODEV);
      rgpulse(pw, t7, 3.0e-6, 5.0e-6);	/* H1 purge pulse */


   status(D);
      txphase(t4);
      genqdphase(t3, xrfdev);
      delay(2*tauxh - pw - 8.0e-6 - PRG_STOP_DELAY - POWER_DELAY);

      gensim2pulse(pw, pwx, t4, t3, 0.0, 5.0e-6, TODEV, xrfdev);
      txphase(zero);
      genqdphase(zero, xrfdev);
      delay(deltaxh - (2*pw/M_PI) - (pwx - pw)/2 - 5.0e-6);

      gensim2pulse(2*pw, 2*pwx, zero, zero, 0.0, 0.0, TODEV, xrfdev);
      delay(deltaxh - rof2 - POWER_DELAY - PRG_START_DELAY);

      rlpower( ((xrfdev == DODEV) ? dpwr : dpwr2), xrfdev);
      delay(rof2);
      rcvron();


   status(E);
      setreceiver(t6);
}
