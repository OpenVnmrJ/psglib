/* refhsq.c         heteronuclear Overbodenhausen experiment using
                    refocused INEPT and broadband decoupling of H1
                    during t1; this sequence has RF duality.


   Parameters:

      sspul = 'y':  selects for Trim(x)-Trim(y) sequence at the start
		    of the pulse sequence
              'n':  normal experiment
        fad = 'y':  TPPI axial-peak displacement
              'n':  standard phasecycle
      f1180 = 'y':  the first t1 point is sampled at half the t1 dwell time
              'n':  the first t1 point is sampled at t1 = 0
     satmode = 'y':  presaturation during relaxation period (satdly) with xmtr
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
	      used; the other delay used in the refocused INEPT
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


static int	phs1[4]  = {1,1,3,3},
		phs2[2]  = {0,2},
		phs3[16] = {0,0,0,0,0,0,0,0,2,2,2,2,2,2,2,2},
		phs4[32] = {0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
			    2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2},
		phs5[64] = {0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
			    0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
			    2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,
			    2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2},
		phs6[32] = {0,2,2,0,0,2,2,0,2,0,0,2,2,0,0,2,
			    2,0,0,2,2,0,0,2,0,2,2,0,0,2,2,0},
		phs7[8]  = {0,0,0,0,1,1,1,1};

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
		f1180[MAXSTR],
		hdshape[MAXSTR],
                *dmx;
   int          phase,
		satmove,
		xrfdev,
		t1_counter;
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
   ni = getval("ni");
   sw1 = getval("sw1");
   sstrim = getval("sstrim");
   hd90 = getval("hd90");
   hdres = getval("hdres");
   hdpwr = getval("hdpwr");
   phase = (int) (getval("phase") + 0.5);
   xrfdev = (int) (getval("xrfdev") + 0.5);

   getstr("sspul", sspul);
   getstr("f1180", f1180);
   getstr("fad", fad);
   getstr("satmode", satmode);
   getstr("hdshape", hdshape);	/* H1 t1 decoupling pattern	*/


/* Load phase tables */
   settable(t1, 4, phs1);
   settable(t2, 2, phs2);
   settable(t3, 16, phs3);
   settable(t4, 32, phs4);
   settable(t5, 64, phs5);
   settable(t6, 32, phs6);
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


/* Check for H1 frequency change */
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
      tssub(t2, 1, 4);


/* FAD phase incrementation */
   if (fad[A] == 'y')
   {
      if (ix == 1)
        d2_init = d2;
      t1_counter = (int) ( (d2 - d2_init)*sw1 + 0.5 );
      if (t1_counter % 2)
      {
         tsadd(t2, 2, 4);
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
      if (satmode[A] == 'y')
      {
         if (satmove)
            offset(satfrq, TODEV);

         rlpower(satpwr, TODEV);
         rgpulse(satdly, zero, 4.0e-5, 0.2e-6);
         if (satmove)
            offset(tof, TODEV);

         rlpower(tpwr, TODEV);
      }


   status(B);
/* does nothing now */


   status(C);
      rcvroff();
      rgpulse(pw, zero, 20.0e-6, 0.0);

      genqdphase(zero, xrfdev);
      txphase(zero);
      delay(deltaxh - pwx);
      gensim2pulse(2*pw, 2*pwx, zero, zero, 0.0, 0.0, TODEV, xrfdev);

      txphase(t1);
      genqdphase(t2, xrfdev);
      delay(deltaxh - pwx);
      gensim2pulse(pw, pwx, t1, t2, 0.0, 0.0, TODEV, xrfdev);

      txphase(zero);
      genqdphase(zero, xrfdev);
      delay( tauxh - pwx - (2*pwx/M_PI) );
      gensim2pulse(2*pw, 2*pwx, zero, zero, 0.0, 0.0, TODEV, xrfdev);

      rlpower(hdpwr, TODEV);
      delay(tauxh - pwx - PRG_START_DELAY - POWER_DELAY);

/* Calculate t1 delay */
      t1evol = d2;
      if (f1180[A] == 'y')
         t1evol += 0.5/sw1;

      obsprgon(hdshape, hd90, hdres);
      xmtron();
      delay(t1evol);
      xmtroff();
      obsprgoff();

      rlpower(tpwr, TODEV);
      rgpulse(pw, t7, 3.0e-6, 0.0);


   status(D);
      txphase(zero);
      genqdphase(t5, xrfdev);
      delay(tauxh - pwx - pw - 3.0e-6 - PRG_STOP_DELAY - POWER_DELAY);

      gensim2pulse(2*pw, 2*pwx, zero, t5, 0.0, 0.0, TODEV, xrfdev);

      txphase(t3);
      genqdphase(t4, xrfdev);
      delay( tauxh - pwx - (2*pwx/M_PI) );

      gensim2pulse(pw, pwx, t3, t4, 0.0, 0.0, TODEV, xrfdev);

      txphase(zero);
      genqdphase(zero, xrfdev);
      delay(deltaxh - pwx);
      gensim2pulse(2*pw, 2*pwx, zero, zero, 0.0, 0.0, TODEV, xrfdev);

      rlpower( ((xrfdev == DODEV) ? dpwr : dpwr2), xrfdev);
      delay(rof2);
      rcvron();
      delay(deltaxh - pwx - POWER_DELAY - PRG_START_DELAY - rof2);


   status(E);
      setreceiver(t6);
}
