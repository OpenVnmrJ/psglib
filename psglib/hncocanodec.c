/* hncocanodec.c  -  HN(CO)CA 3D experiment by Kay et al.; uses a
		    frequency-shifted pulse on the Ca spins; has
		    RF duality

                    no proton decoupling is used(as opposed to hncoca.c)

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
     satmode = 'y':  H1 presaturation during relaxtion delay
     satfrq = frequency of 1H presaturation for all periods
     satdly = saturation time during the relaxation period
     satpwr = saturation power for all periods of presaturation with xmtr
         hs = 'yn':  homospoil pulse (hst) during the d1 relaxation delay
       tpwr = power level for 1H transmitter pulses
         pw = 90 degree xmtr pulse length for protons (the observed nucleus)
       dof2 = should be set to C0 frequency
    pwcalvl = power level for Ca decoupler pulses
       pwca = 90 degree decoupler pulse length for Ca at `pwcalvl`
    cashape = frequency-shifted pulse for Ca excitation
    pwcolvl = power level for C0 decoupler pulses
       pwco = 90 degree decoupler pulse length for C0 at `pwcolvl`
   pwco2lvl = power level for shorter CO pulse
      pwco2 = 90 degree decoupler pulse length for C0 at `pwco2lvl`
   pwn15lvl = power level for N15 decoupler pulses
      pwn15 = 90 degree decoupler pulse length for N15
       dpwr = power level for N15 broadband decoupling
       tau1 = first delay in sequence (NH scalar coupling)
       tau2 = second delay in sequence
       tau3 = third delay in sequence
       tau4 = fourth delay in sequence
         dm = 'nnnnn':  no broadband decoupling of N15 during acquisition
              'nnnny':  broadband heteronuclear decoupling of N15
		        during acquisition
        dm2 = 'n':  no broadband decoupling of C13 during acquisition

      phase = 1,2:  hypercomplex experiment with F1 quadrature (complex F1-FT)
     phase2 = 1,2:  hypercomplex experiment with F2 quadrature (complex F2-FT)
                modified from hncoca1FSP.c(S.Farmer,Varian)  1/27/92
*/



#include <standard.h>
#include <math.h>

#define MIN_DELAY	0.2e-6		/* shortest executable delay */


static int	 phs1[4] = {0,0,2,2},
		 phs2[2] = {0,2},
		 phs3[8] = {0,0,0,0,2,2,2,2},
		phs5[16] = {0,0,0,0,0,0,0,0,2,2,2,2,2,2,2,2},
		 rec[16] = {0,2,2,0,2,0,0,2,2,0,0,2,0,2,2,0};

static double	d2_init = 0.0,
		d3_init = 0.0;


pulsesequence()
{
/* VARIABLE DECLARATION */
   char         satmode[MAXSTR],
		sspul[MAXSTR],
		fad1[MAXSTR],
		fad2[MAXSTR],
		f1180[MAXSTR],
		f2180[MAXSTR],
		cashape[MAXSTR];
   int          phase,
		phase2,
		satmove,
		t1_counter,
		t2_counter,
		c13dev = DODEV,
		n15dev = DO2DEV;
   double       ss,
		t1evol,
		t2evol_1,
		t2evol_2,
		sw1,
		sw2,
                pwcalvl,
                pwca,
                pwcolvl,
                pwco,
                pwco2lvl,
                pwco2,
                pwn15lvl,
                pwn15,
                satfrq,
                satdly,
                satpwr,
		tau1,
		tau2,
		tau3,
		tau4;


/* Load variables */
   satfrq = getval("satfrq");
   satdly = getval("satdly");
   satpwr = getval("satpwr");
   pwcalvl = getval("pwcalvl");
   pwca = getval("pwca");
   pwcolvl = getval("pwcolvl");
   pwco = getval("pwco");
   pwco2lvl = getval("pwco2lvl");
   pwco2 = getval("pwco2");
   pwn15lvl = getval("pwn15lvl");
   pwn15 = getval("pwn15");
   tau1 = getval("tau1");
   tau2 = getval("tau2");
   tau3 = getval("tau3");
   tau4 = getval("tau4");
   ss = getval("ss");
   sw1 = getval("sw1");
   sw2 = getval("sw2");
   phase = (int) (getval("phase") + 0.5);
   phase2 = (int) (getval("phase2") + 0.5);

   getstr("cashape", cashape);
   getstr("sspul", sspul);
   getstr("fad1", fad1);
   getstr("fad2", fad2);
   getstr("satmode", satmode);
   getstr("f1180", f1180);
   getstr("f2180", f2180);


/* Load phase cycles */
   settable(t1, 4, phs1);
   settable(t2, 2, phs2);
   settable(t3, 8, phs3);
   settable(t4, 16, rec);
   settable(t5, 16, phs5);


/* Check conditions */
   satmove = ( fabs(tof - satfrq) >= 0.1 );

   if ( (dm2[A] == 'y') || (dm2[B] == 'y') || (dm2[C] == 'y')
		|| (dm2[D] == 'y') )
   {
      text_error("`dm2` must be to 'nnnnn' or 'nnnny' for N15 decoupling\n");
      abort(1);
   }

   if ( (dm[A] == 'y') || (dm[B] == 'y') || (dm[C] == 'y')
		|| (dm[D] == 'y') || (dm[E] == 'y') )
   {
      text_error("`dm` must be to 'n' for C13\n");
      abort(1);
   }


/* Determine steady-state mode */
   if (ss < 0)
   {
      ss *= (-1);
      initval(ss, ssval);
      initval(ss, ssctr);
   }


/* Add in States-Haberkorn element */
   if (phase2 == 2)		/* N15 t1 element */
      tsadd(t2, 1, 4);

   if (phase == 2)		/* C13 t2 element */
      tsadd(t3, 1, 4);


/* Add in FAD */
   if (fad2[A] == 'y')		/* for N15 */
   {
      if (ix == 1)
        d3_init = d3;

      t2_counter = (int) ( (d3 - d3_init)*sw2 + 0.5 );
      if (t2_counter % 2)
      {
         tsadd(t2, 2, 4);	/* first N15 90-degree pulse 	*/
         tsadd(t4, 2, 4);	/* receiver phase cycle		*/
      }
   }

   if (fad1[A] == 'y')		/* for C13 */
   {
      if (ix == 1)
        d2_init = d2;
 
      t1_counter = (int) ( (d2 - d2_init)*sw1 + 0.5 );
      if (t1_counter % 2)
      {   
         tsadd(t3, 2, 4);       /* first C13a 90-degree pulse   */
         tsadd(t4, 2, 4);       /* receiver phase cycle         */
      }
   }


/* BEGIN ACTUAL PULSE SEQUENCE CODE */
   status(A);
      rlpower(tpwr, TODEV);		/*  H1 hard-pulse power level	*/
      rlpower(pwn15lvl, n15dev);	/* N15 hard-pulse power level	*/
      rlpower(pwcolvl, c13dev); 	/* CO hard-pulse power		*/

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
   status(B);
      rcvroff();

      rgpulse(pw, zero, rof1, 0.0);
      genqdphase(t1, n15dev);		/* N15 phase */
      genqdphase(zero, c13dev);		/* C13 phase */
      delay(tau1);

      gensim2pulse(2*pw, 2*pwn15, zero, t1, 0.0, 0.0, TODEV, n15dev);
      txphase(t1);			/* H1 phase  */
      genqdphase(t2, n15dev);		/* N15 phase */
      delay(tau1);

      getelem(t1, ct, v14);
      add(one, v14, v14);

      gensim2pulse(pw, pwn15, v14, t2, 0.0, 0.0, TODEV, n15dev);
      txphase(zero);
      genqdphase(zero, n15dev);
      delay(tau2);

      t1evol = d3;	/* N15 is governed by ni2 and sw2 */
      if (f2180[0] == 'y')
         t1evol += 0.5/sw2;

      if (t1evol > MIN_DELAY)
      {
         t1evol -= 2*pw;
         if (t1evol < MIN_DELAY)
            t1evol = 0.0;
      }

      delay(t1evol/2);
      rgpulse(2*pw, zero, 0.0, 0.0);
      delay(t1evol/2);


   status(C);
      delay(tau3);
      gensim2pulse(2*pwn15, 2*pwco, zero, zero, 0.0, 0.0, n15dev, c13dev);
      delay(tau2 + tau3);
      gensim2pulse(pwn15, pwco, zero, zero, 0.0, 0.0, n15dev, c13dev);

      delay(tau4/2);
      rlpower(pwcalvl, c13dev);
      genqdphase(t3, c13dev);
      delay( (tau4/2) - POWER_DELAY - WFG_START_DELAY );
      genshaped_pulse(cashape, pwca, t3, 0.0, 0.0, 0.0, 0.0, c13dev);

      t2evol_1 = t2evol_2 = d2/2;	/* C13 is governed by ni and sw */
      if (f1180[0] == 'y')
      {
         t2evol_1 += 0.25/sw1;
         t2evol_2 += 0.25/sw1;
      }

      if ( (t2evol_1 > MIN_DELAY) || (t2evol_2 > MIN_DELAY) )
      {
         t2evol_1 -= (2*pwca/M_PI) + pwco2 + POWER_DELAY + WFG_STOP_DELAY;
         t2evol_2 -= (2*pwca/M_PI) + pwco2 + POWER_DELAY + WFG_START_DELAY;

         if (t2evol_1 < MIN_DELAY)
            t2evol_1 = 0.0;
         if (t2evol_2 < MIN_DELAY)
            t2evol_2 = 0.0;
      }

      genqdphase(zero, c13dev);
      rlpower(pwco2lvl, c13dev);
      delay(t2evol_1);

      gensim2pulse(2*pw, 2*pwco2, zero, zero, 0.0, 0.0, TODEV, c13dev);

      rlpower(pwcalvl, c13dev);
      delay(t2evol_2);

      genshaped_pulse(cashape, pwca, zero, 0.0, 0.0, 0.0, 0.0, c13dev);
      delay( (tau4/2) - WFG_STOP_DELAY );
      rlpower(pwcolvl, c13dev);
      delay( (tau4/2) - POWER_DELAY );


   status(D);
      gensim2pulse(pwn15, pwco, t5, zero, 0.0, 0.0, n15dev, c13dev);
      delay(tau3 + tau2);
      gensim2pulse(2*pwn15, 2*pwco, zero, zero, 0.0, 0.0, n15dev, c13dev);
      delay(tau3);

      rgpulse(2*pw, zero, 0.0, 0.0);
      delay(tau2);
      gensim2pulse(pw, pwn15, zero, zero, 0.0, 0.0, TODEV, n15dev);
      delay(tau1);
      gensim2pulse(2*pw, 2*pwn15, zero, zero, 0.0, 0.0, TODEV, n15dev);

      rlpower( ((n15dev == DODEV) ? dpwr : dpwr2), n15dev);
      rlpower( ((c13dev == DODEV) ? dpwr : dpwr2), c13dev);
      delay(rof2);
      rcvron();
      delay(tau1 - rof2 - 2*POWER_DELAY);


   status(E);
      setreceiver(t4);
}
