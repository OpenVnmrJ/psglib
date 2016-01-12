/*
 hncoca          -  HN(CO)CA 3D experiment by Kay ; uses a
	   	    frequency-shifted pulse on the Ca spins
                    to permit use of three channels instead of four.

         ref: A.Bax and M.Ikura, J. Biomol.NMR,1,99(1991)


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
       dof  = should be set to C0 frequency
      dof2  = should be set in center of 15N NH region
       pwc13 = 90 degree decoupler pulse length for C13 at `pwc13lvl`
    cashape = frequency-shifted pulse for Ca excitation

               create cashape entry using "convolute" macro, for example
                convolute('275hard',ca90,55,-18000) for a 55us 90-degree
                pulse on the alpha carbons, 18.3 kHz upfield from the C13
                xmtr position. In this case, "275hard" is a shapelib file
                containing 275 identical lines of "0.0  1023.0"
                (275x0.2us)=55us

    pwc13lvl = power level for Ca decoupler pulses
   pwn15lvl = power level for N15 decoupler pulses
      pwn15 = 90 degree decoupler pulse length for N15
       dpwr = power level for N15 broadband decoupling
       tau1 = first delay in sequence (NH scalar coupling)
       tau2 = second delay in sequence
       tau3 = third delay in sequence
       tau4 = fourth delay in sequence
        dm2 = 'nnn':  no broadband decoupling of N15 during acquisition
              'nny':  broadband heteronuclear decoupling of N15
        dm  = 'n':  no broadband decoupling of C13 during acquisition
     lsfrq1 = same amount as used to create cashape file
      phase = 1,2:  hypercomplex experiment with F1 quadrature (complex F1-FT)
     phase2 = 1,2:  hypercomplex experiment with F2 quadrature (complex F2-FT)



*/

#include <standard.h>
#include <math.h>

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
		t2_counter;
   double       ss,
		t1evol,
		t2evol_1,
		t2evol_2,
		sw1,
		sw2,
                pwc13lvl,
                pwc13,
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
   pwc13 = getval("pwc13");
   pwc13lvl = getval("pwc13lvl");
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

   if ( (dm2[A] == 'y') || (dm2[B] == 'y') )
   {
      text_error("`dm2` must be to 'nnn' or 'nny' for N15 decoupling\n");
      abort(1);
   }

   if ( (dm[A] == 'y') || (dm[B] == 'y') || (dm[C] == 'y'))
   {
      text_error("`dm` must be to 'nnn' for C13\n");
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
      rlpower(pwn15lvl,DO2DEV);	        /* N15 hard-pulse power level	*/
      rlpower(pwc13lvl, DODEV);    	/* C13 hard-pulse power		*/

      if (sspul[A] == 'y')
      {
         rgpulse(200*pw, zero, rof1,0.0);
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
      delay(tau1);

      sim3pulse(2*pw,0.0, 2*pwn15, zero,zero, t1, 0.0, 0.0);
      delay(tau1);

      getelem(t1, ct, v14);
      add(one, v14, v14);

      sim3pulse(pw, 0.0, pwn15, v14,zero, t2, 0.0, 0.0);
      delay(tau2);

      t1evol = d3;	/* N15 is governed by ni2 and sw2 */
      if (f2180[0] == 'y')
         t1evol += 0.5/sw2;

      delay(t1evol/2);
      rgpulse(2*pw, zero, 0.0, 0.0);
      delay(t1evol/2);
      rlpower(pwc13lvl+6,DODEV);
      delay(tau3);
      sim3pulse(0.0, pwc13, 2*pwn15,zero, zero, zero, 0.0, 0.0);
      rlpower(pwc13lvl,DODEV);
      delay(tau2 + tau3);
      sim3pulse(0.0,pwc13,pwn15,zero, zero, zero, 0.0, 0.0);
      delay(tau4);
      decshaped_pulse(cashape, pwc13, t3, 0.0, 0.0);

      t2evol_1 = t2evol_2 = d2/2;	/* C13 is governed by ni and sw */
      if (f1180[0] == 'y')
      {
         t2evol_1 += 0.25/sw1;
         t2evol_2 += 0.25/sw1;
      }

      rlpower(pwc13lvl+6, DODEV);
      delay(t2evol_1);
      simpulse(2*pw, pwc13, zero, zero, 0.0, 0.0);
      rlpower(pwc13lvl, DODEV);
      delay(t2evol_2);
      decshaped_pulse(cashape, pwc13, zero, 0.0, 0.0);
      delay( (tau4/2) );
      rlpower(pwc13lvl, DODEV);
      delay( (tau4/2) );
      sim3pulse(0.0,pwc13,pwn15,zero, zero,t5, 0.0, 0.0);
      rlpower(pwc13lvl+6,DODEV);
      delay(tau3 + tau2);
      sim3pulse(0.0,pwc13,2*pwn15,zero, zero, zero, 0.0, 0.0);
      rlpower(pwc13lvl,DODEV);
      delay(tau3);
      rgpulse(2*pw, zero, 0.0, 0.0);
      delay(tau2);
      sim3pulse(pw,0.0, pwn15, zero,zero, zero, 0.0, 0.0);
      delay(tau1);
      sim3pulse(2*pw,0.0, 2*pwn15,zero, zero, zero, 0.0, 0.0);

      rlpower( dpwr,DODEV);
      rlpower( dpwr2,DO2DEV);
      delay(rof2);
      rcvron();
      delay(tau1 );


   status(C);
      setreceiver(t4);
}
