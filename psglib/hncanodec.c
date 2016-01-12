/* hncanodec.c  -  modified HNCA experiment by Kay et al.; uses shifted
		  laminar pulses to selectively excite the Ca spins
		  while on resonance with the CO spins; on-resonance
		  selective decoupling of the CO spins is employed during
		  tau(cn) and t1; refocusing of the Ca and CO spins is
		  used during t2; this sequence has RF duality.

                       no proton decoupling is used(as in hnca.c)


   Parameters:

      sspul = 'y':  selects for Trim(x)-Trim(y) sequence at the start
		    of the pulse sequence
      f1180 = 'n':  standard t1 timing
              'y':  modified t1 timing for t1(1) = half the dwell time
      f2180 = 'n':  standard t2 timing
              'y':  modified t2 timing for t2(1) = half the dwell time
              'n':  normal experiment
       fad1 = 'y':  TPPI axial-peak displacement along t1
              'n':  standard phasecycle
       fad2 = 'y':  TPPI axial-peak displacement along t2 (3D experiment)
              'n':  standard phasecycle
     satmode = 'nnnn':  no H1 presaturation
              'ynnn':  H1 presaturation during relaxation delay
              'nynn':  H1 presaturation during the first tau(cn) period
              'nnyn':  H1 presaturation during the t1 period
              'nnny':  H1 presaturation during the second tau(cn) period
     satfrq = frequency of H1 presaturation for all periods
     satdly = saturation time during the relaxation period
     satpwr = saturation power for all periods of presaturation with xmtr
         hs = 'yn':  homospoil pulse (hst) during the d1 relaxation delay
       tpwr = power level for H1 transmitter pulses
         pw = 90 degree xmtr pulse length for protons (the observed nucleus)
    cashape = pattern for shifted laminar pulse on the Ca spins
    pwcalvl = power level for Ca decoupler pulses
       pwca = 90 degree decoupler pulse length for Ca at `pwcalvl`
      cafrq = frequency for Ca spins
    pwcolvl = power level for selective C0 decoupling
       pwco = 90 degree pulse length for selective C0 decoupling at `pwcolvl`
    coshape = decoupling pattern for selective CO decoupling
      cores = tip-angle resolution for selective CO decoupling
      cofrq = frequency for CO spins
     c13dev = RF device for C13
   pwn15lvl = power level for N15 decoupler pulses
      pwn15 = 90 degree decoupler pulse length for N15
     n15dev = RF device for N15
        jnh = one-bond heteronuclear coupling constant to NH (in Hz)
        jcn = one-bond heteronuclear coupling constant to CN (in Hz)
  dpwr(N15) = power level for N15 broadband decoupling
  dpwr(C13) = power level for C13 broadband decoupling
    dm(N15) = 'nn':  no broadband decoupling of N15 during acquisition
              'ny':  broadband heteronuclear decoupling of N15 during
		     acquisition
    dm(C13) = 'nn':  no broadband decoupling of C13 during acquisition

      phase = 1,2:  hypercomplex experiment with F1 quadrature (complex F1-FT)
     phase2 = 1,2:  hypercomplex experiment with F2 quadrature (complex F2-FT)
                  modified from hnca1FSP.c(S.Farmer,Varian)  1/27/92
*/



#include <standard.h>
#include <math.h>

#define MIN_DELAY	0.2e-6		/* shortest executable delay */
#define MIN_J		0.1		/* Hz  */
#define MIN_NULL	0.0001		/* sec */


static int	   count = 0,
		phs1[4]  = {0,0,2,2},
		phs2[2]  = {0,2},
		phs3[2]  = {0,1},
		phs4[2]  = {0,2},
		phs5[2]  = {0,2},
		phs6[2]  = {0,1},
		phs7[64] = {0,2,2,0,2,0,0,2,2,0,0,2,0,2,2,0,
			    2,0,0,2,0,2,2,0,0,2,2,0,2,0,0,2,
			    2,0,0,2,0,2,2,0,0,2,2,0,2,0,0,2,
			    0,2,2,0,2,0,0,2,2,0,0,2,0,2,2,0};

static double	d2_init = 0.0,
		d3_init = 0.0,
		t1adj_1 = 0.0,
		t1adj_2 = 0.0,
		tmpd2 = 0.0;

extern double	getvalnwarn();


/*-------------------------------
|				|
|       pulsesequence()/0	|
|				|
+------------------------------*/
pulsesequence()
{
/* VARIABLE DECLARATION */
   char         satmode[MAXSTR],
		sspul[MAXSTR],
		fad1[MAXSTR],
		f1180[MAXSTR],
		f2180[MAXSTR],
		prtinfo[MAXSTR],
		cashape[MAXSTR],
		coshape[MAXSTR],
		*dmc13,
		*dmn15;
   int          phase,
		satmove,
		t1_counter,
		exp3D,
		c13dev,
		n15dev;
   double       ss,
		sw1,
		sw2,
		delayval,
		t1evol_1,
		t1evol_2,
		d2local,
		cafrq,
		cofrq,
                pwcalvl,
                pwca,
                pwcolvl,
                pwco,
		cores,
                pwn15lvl,
                pwn15,
                jnh,
                jcn,
		deltanh,
                deltacn,
                satfrq,
                satdly,
                satpwr;


/* LOAD VARIABLES */
   satfrq = getval("satfrq");
   satdly = getval("satdly");
   satpwr = getval("satpwr");
   cofrq = getval("cofrq");
   cafrq = getval("cafrq");
   pwcalvl = getval("pwcalvl");
   pwca = getval("pwca");
   pwcolvl = getval("pwcolvl");
   pwco = getval("pwco");
   cores = getval("cores");
   pwn15lvl = getval("pwn15lvl");
   pwn15 = getval("pwn15");
   jnh = getval("jnh");
   jcn = getval("jcn");
   ss = getval("ss");
   sw1 = getval("sw1");
   phase = (int) (getval("phase") + 0.5);
   c13dev = (int) (getval("c13dev") + 0.5);
   n15dev = (int) (getval("n15dev") + 0.5);

   getstr("cashape", cashape);
   getstr("coshape", coshape);
   getstr("sspul", sspul);
   getstr("fad1", fad1);
   getstr("satmode", satmode);
   getstr("f1180", f1180);
   getstr("prtinfo", prtinfo);


/* INITIALIZE VARIABLES */
   deltanh = ( (jnh > MIN_J) ? 1/(2*jnh) : 0.0 );
   deltacn = ( (jcn > MIN_J) ? 1/(2*jcn) : 0.0 );


/* LOAD PHASE TABLE */
   settable(t1, 4, phs1);
   settable(t2, 2, phs2);
   settable(t3, 2, phs3);
   settable(t4, 2, phs4);
   settable(t5, 2, phs5);
   settable(t6, 2, phs6);
   settable(t7, 64, phs7);

   setdivnfactor(t3, 32);
   setdivnfactor(t4, 4);
   setdivnfactor(t5, 16);
   setdivnfactor(t6, 8);


/* CHECK CONDITIONS */
   satmove = ( fabs(tof - satfrq) >= 0.1 );
   exp3D = ( (int) (getvalnwarn("ni2") + 0.5) > 0 );

   if ( (n15dev != DODEV) && (n15dev != DO2DEV) )
   {
      text_error("invalid RF device for N15\n");
      abort(1);
   }
   else if ( (c13dev != DODEV) && (c13dev != DO2DEV) ) 
   { 
      text_error("invalid RF device for C13\n");
      abort(1); 
   }
   else if (c13dev == n15dev)
   {
      text_error("N15 and C13 RF devices must be different\n");
      abort(1);
   }

   if (satmode[C] == 'y')
   {
      text_error("`satmode` must be 'n' during status C\n");
      abort(1);
   }

   dmn15 = ( (n15dev == DODEV) ? dm : dm2 );
   dmc13 = ( (c13dev == DODEV) ? dm : dm2 );

   if ( (strcmp(dmn15, "n") != 0) && (strcmp(dmn15, "nn") != 0) &&
	(strcmp(dmn15, "ny") != 0) )
   {
      text_error("`dm(N15)` must be to 'n', 'nn', or 'ny'\n");
      abort(1);
   }

   if ( (strcmp(dmc13, "n") != 0) && (strcmp(dmc13, "nn") != 0) )
   {
      text_error("`dm(C13)` must be to 'n' or 'nn'\n");
      abort(1);
   }


/* DETERMINE STEADY-STATE MODE */
   if (ss < 0)
   {
      ss *= (-1);
      initval(ss, ssval);
      initval(ss, ssctr);
   }


/* ADD IN STATES-HABERKORN ELEMENT */
   if (phase == 2)
      tsadd(t2, 1, 4);

   if (exp3D)
   {
      int	phase2;

      phase2 = (int) (getval("phase2") + 0.5);
      if (phase2 == 2)
         tssub(t1, 1, 4);
   }


/* ADD IN FAD */
   if (fad1[A] == 'y')
   {
      if (ix == 1)
        d2_init = d2;

      t1_counter = (int) ( (d2 - d2_init)*sw1 + 0.5 );
      if (t1_counter % 2)
      {
         tsadd(t2, 2, 4);	/* first C13 90 phase cycle	*/
         tsadd(t7, 2, 4);	/* receiver phase cycle		*/
      }
   }

   if (exp3D)
   {
      char	fad2[MAXSTR];
      int	t2_counter;

      sw2 = getval("sw2");
      getstr("fad2", fad2);
      getstr("f2180", f2180);

      if (fad2[A] == 'y')
      {
         if (ix == 1)
           d3_init = d3;
 
         t2_counter = (int) ( (d3 - d3_init)*sw2 + 0.5 );
         if (t2_counter % 2)
         {   
            tsadd(t1, 2, 4);      /* first N15 90 phase cycle     */
            tsadd(t7, 2, 4);      /* receiver phase cycle         */
         }
      }
   }


/* BEGIN ACTUAL PULSE SEQUENCE CODE */
   status(A);
      rlpower(tpwr, TODEV);		/*  H1 hard-pulse power level	*/
      rlpower(pwn15lvl, n15dev);	/* N15 hard-pulse power level	*/
      rlpower(pwcolvl, c13dev); 	/* CO pulse power		*/
      offset(cofrq, c13dev); 		/* CO frequency			*/

      if (sspul[A] == 'y')
      {
         rgpulse(200*pw, zero, rof1, 0.0);
         rgpulse(200*pw, one, 0.0, rof2);
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


/* Pulse train */
      rcvroff();
      genqdphase(zero, c13dev);         /* CO selective decoupling phase */
      (void) prg_dec_on(coshape, pwco, cores, c13dev);
					/* start PRG-dec pattern running */

/* INEPT transfer from H1 to N15 */
      rgpulse(pw, zero, rof1, 5.0e-6);
      txphase(zero);			/* H1 phase  */
      genqdphase(t1, n15dev);		/* N15 phase */

      delay(deltanh - 5.0e-6);

      gensim2pulse(2*pw, pwn15, zero, t1, 0.0, 5.0e-6, TODEV, n15dev);
      rfon(c13dev);			/* turn on CO RF	*/
      txphase(one);			/* H1 phase  		*/
      genqdphase(t3, n15dev);		/* N15 phase 		*/
      delay(deltanh - 5.0e-6);

      rgpulse(pw, one, 0.0, 0.0);


/* HMQC established between N15 and C13 */
      if (satmode[B] == 'y')
      {
         rlpower(satpwr, TODEV);
         txphase(zero);			/* H1 phase for presaturation */

         if (satmove)
         {
            offset(satfrq, TODEV);
            rgpulse(deltacn - 2*(POWER_DELAY + OFFSET_DELAY) - POWER_DELAY
			- WFG_START_DELAY - PRG_STOP_DELAY - deltanh - pw,
			zero, 0.0, 0.0);
            offset(tof, TODEV);
         }
         else
         {
            rgpulse(deltacn - 3*POWER_DELAY - WFG_START_DELAY - PRG_STOP_DELAY
			 - deltanh - pw, zero, 0.0, 0.0);
         }

         rlpower(tpwr, TODEV);
      }
      else
      {
         delay(deltacn - deltanh - pw - POWER_DELAY - WFG_START_DELAY
			- PRG_STOP_DELAY);
         txphase(zero);
      }


/* Calculate t1 evolution times */
      if (ix == 1)
      {
         t1adj_1 = (2*pwca/M_PI) + pw + WFG_START_DELAY + PRG_STOP_DELAY
			+ POWER_DELAY;
         t1adj_2 = pwn15 + pw;

         t1evol_1 = (8 * sw1) * ( MIN_DELAY + t1adj_1 );
         t1evol_2 = (8 * sw1) * ( MIN_DELAY + t1adj_2 );
         tmpd2 = ( (t1evol_1 > t1evol_2) ? t1evol_1 : t1evol_2 );

         if (prtinfo[0] == 'y')
         {
            (void) printf( "\nMaximum C13 spectral width = %lf\n",
		( (f1180[0] == 'y') ? (sw1/tmpd2) : (2*sw1)/tmpd2 ) );
         }

         count = (int) (tmpd2) + 1;
         tmpd2 = count * (0.5/sw1);
         if (prtinfo[0] == 'y')
         {
            char msge[128];

            (void) sprintf(msge, "count = %d\n", count);
            text_error(msge);
         }
      }


      if ( (f1180[0] != 'y') && (count < 3) )
      {
         t1evol_1 = (d2/4) - t1adj_1;
         t1evol_2 = (d2/4) - t1adj_2;
      }
      else
      {
         if ( (count > 2) || ((f1180[0] == 'y') && (count > 1)) )
         {
            if (ix == 1)
               text_error("WARNING:  sw1 is too large for proper t1 timing\n");
         }

         d2local = d2 + tmpd2;
         t1evol_1 = (d2local/4) - t1adj_1;
         t1evol_2 = (d2local/4) - t1adj_2;
      }

      if (t1evol_1 < MIN_DELAY)
         t1evol_1 = 0.0;

      if (t1evol_2 < MIN_DELAY)
         t1evol_2 = 0.0;


/* Start C13 pulsing */
      rfoff(c13dev);			/* turn off CO RF		*/
      rlpower(pwcalvl, c13dev);		/* Ca pulse power level		*/
      genqdphase(t2, c13dev);		/* Ca phase			*/
      prg_dec_off(2, c13dev);		/* shutdown C13 WFG		*/

      genshaped_pulse(cashape, pwca, t2, 0.0, 0.0, 0.0, 0.0, c13dev);

      rlpower(pwcolvl, c13dev);		/* CO decoupling power level	*/
      genqdphase(zero, c13dev);		/* CO decoupling phase		*/
      (void) prg_dec_on(coshape, pwco, cores, c13dev);
					/* turn on CO decoupling	*/
      rfon(c13dev);			/* turn on CO RF		*/

      delay(t1evol_1);
      rgpulse(2*pw, zero, 0.0, 0.0);
      delay(t1evol_2);

      genrgpulse(2*pwn15, t3, 0.0, 0.0, n15dev);

      delay(t1evol_2);
      rgpulse(2*pw, zero, 0.0, 0.0);
      delay(t1evol_1);

      rfoff(c13dev);			/* turn off CO RF		*/
      rlpower(pwcalvl, c13dev);      	/* Ca pulse power level  	*/
      genqdphase(t4, c13dev);		/* Ca phase			*/
      genqdphase(zero, n15dev);		/* N15 phase			*/
      prg_dec_off(2, c13dev);		/* shutdown C13 WFG		*/

      genshaped_pulse(cashape, pwca, t4, 0.0, 0.0, 0.0, 0.0, c13dev);

      rlpower(pwcolvl, c13dev);		/* CO decoupling power level	*/
      genqdphase(zero, c13dev);		/* CO decoupling phase		*/
      (void) prg_dec_on(coshape, pwco, cores, c13dev);
      rfon(c13dev);			/* turn on CO RF		*/


/* Remove HMQC between N15 and C13 */
      delayval = deltacn - deltanh - POWER_DELAY - WFG_STOP_DELAY -
		   PRG_STOP_DELAY;
      if (exp3D)
         delayval -= 8*pwca + OFFSET_DELAY + POWER_DELAY + 5.0e-6;

      if (satmode[D] == 'y')
      {
         txphase(zero);			/* H1 phase			*/
         rlpower(satpwr, TODEV);	/* H1 presaturation power level	*/

         if (satmove)
         {
            offset(satfrq, TODEV);
            rgpulse(delayval - 2*POWER_DELAY - 2*OFFSET_DELAY, zero, 0.0, 0.0);
            offset(tof, TODEV);
         }
         else
         {
            rgpulse(delayval - 2*POWER_DELAY, zero, 0.0, 0.0);
         }

         ( (exp3D) ? txphase(zero) : txphase(t5) );
         rlpower(tpwr, TODEV);
      }
      else
      {
         ( (exp3D) ? txphase(zero) : txphase(t5) );
         delay(delayval);
      }


/* Turn off all CO decoupling */
      rfoff(c13dev);			/* turn off CO RF		*/
      prg_dec_off(2, c13dev);		/* shutdown C13 WFG		*/


/* 3D section for N15 t2 evolution */
      if (exp3D)
      {
         double	t2evol;

         t2evol = d3;
         if (f2180[0] == 'y')
            t2evol += (0.5/sw2);

         offset(cafrq, c13dev);		/* Ca center frequency		*/
         rlpower(pwcalvl, c13dev);	/* C13 broadband pulse power 	*/
         delay(5.0e-6);

         delay(t2evol/2);

         genrgpulse(pwca, zero, 0.0, 0.0, c13dev);
         gensim2pulse(2*pw, 2*pwca, zero, one, 0.0, 0.0, TODEV, c13dev);
         genrgpulse(pwca, zero, 0.0, 0.0, c13dev);

         txphase(t5);			/* H1 phase	*/
         delay(t2evol/2);

         genrgpulse(pwca, zero, 0.0, 0.0, c13dev);
         genrgpulse(2*pwca, one, 0.0, 0.0, c13dev);
         genrgpulse(pwca, zero, 0.0, 0.0, c13dev);
      }


/* INEPT transfer from N15 back to H1 */
      rgpulse(pw, t5, 0.0, 5.0e-6);
      txphase(t6);			/* H1 phase  */
      delay(deltanh - pw - 5.0e-6);

      gensim2pulse(2*pw, pwn15, t6, zero, 0.0, 0.0, TODEV, n15dev);

      rlpower(dpwr, DODEV);		/* N15 or C13 decoupling power */
      rlpower(dpwr2, DO2DEV);		/* C13 or N15 decoupling power */
      genqdphase(zero, n15dev);		/* N15 decoupling phase */
      genqdphase(zero, c13dev);		/* C13 final phase */
      delay(rof2);
      rcvron();

      delay(deltanh - 2*POWER_DELAY - PRG_START_DELAY - rof2);


/* Start N15 broadband decoupling */
   status(B);
      setreceiver(t7);
}
