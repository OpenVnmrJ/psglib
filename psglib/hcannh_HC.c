/* hcannh_HC  -

	relayed Ha-->Ca-->N-->Hn 4D experiment with CT t2 and t3
	evolution times; uses SLP pulses to selectively excite the
	CO resonances; the N15 RF and C13 RF channels are user-
	selectable.

	Ha = t1 dimension	( ni, sw1)
	Ca = t2 dimension	(ni2, sw2)
	 N = t3 dimension	(ni3, sw3)
	Hn = t4 dimension	( np, sw )


   Parameters:

      sspul = 'y':  selects for Trim(x)-Trim(y) sequence at the start
		    of the pulse sequence
     sstrim = length of Trim pulse for sspul = 'y'
      f1180 = 'n':  standard t1 timing
              'y':  modified t1 timing for t1(1) = half the dwell time
      f2180 = 'n':  standard t2 timing
              'y':  modified t2 timing for t2(1) = half the dwell time
      f3180 = 'n':  standard t3 timing
              'y':  modified t3 timing for t3(1) = half the dwell time
       fad1 = 'y':  TPPI axial-peak displacement along t1
              'n':  standard phasecycle
       fad2 = 'y':  TPPI axial-peak displacement along t2 (3D experiment)
              'n':  standard phasecycle
       fad3 = 'y':  TPPI axial-peak displacement along t3 (4D experiment)
              'n':  standard phasecycle
      scuba = length of delay used in SCUBA pre-sequence; if scuba > 200 ns,
	      the SCUBA pre-sequence is performed
     satflg = 'y':  H1 presaturation during relaxtion delay
     satfrq = frequency of 1H presaturation for all periods
     satdly = saturation time during the relaxation period
     satpwr = saturation power for all periods of presaturation with xmtr
      hdlvl = power level for H1 broadband decoupling during tau(cn),
	      t1, and t2
    hdshape = H1 broadband decoupling sequence
      hdres = tip-angle resolution for H1 broadband decoupling sequence
       hd90 = 90 degree pulse length for H1 broadband decoupling
      hdfrq = NOT DEFINED ==> tof is assumed to define the center of the H1
              broadband decoupling profile
         hs = 'yn':  homospoil pulse (hst) during the d1 relaxation delay
        jch = one-bond H-Ca heteronuclear coupling constant (in Hz)
    deltach = if deltach < 0.2e-6, deltach = 1/(2*jch)
      tauch = should be set to 3.57 ms (??) (or 1/[2*140]); if tauch < 0.2e-6,
              taunh = 1/(2*jch)
        jnh = one-bond H-N heteronuclear coupling constant (in Hz)
    deltanh = if deltanh < 0.2e-6, deltanh = 1/(2*jnh)
      taunh = should be set to 5.55 ms (or 1/[2*90]); if taunh < 0.2e-6,
              taunh = 1/(2*jnh)
        jcn = effective heteronuclear N-Ca coupling constant (in Hz)
       tpwr = power level for 1H transmitter pulses
         pw = 90 degree xmtr pulse length for protons (the observed nucleus)
      cafrq = frequency for Ca spins (not explicitly used in pulse sequence)
      cofrq = frequency for CO spins
   pwc13lvl = power level for hard C13 decoupler pulses
      pwc13 = 90 degree decoupler pulse length for C13 at `pwc13lvl`
    pwcalvl = power level for Ca decoupler pulses
       pwca = 90 degree decoupler pulse length for Ca at `pwcalvl`
    coshape = frequency-shifted pulse for CO excitation
    pwcolvl = power level for 90-degree C0 decoupler pulses
       pwco = 90 degree decoupler pulse length for C0 at `pwcolvl`
   pwn15lvl = power level for N15 decoupler pulses
      pwn15 = 90 degree decoupler pulse length for N15
  dpwr(N15) = power level for N15 broadband decoupling
  dpwr(C13) = power level for C13 broadband decoupling
    dm(N15) = 'nn':  no broadband decoupling of N15 during acquisition
              'ny':  broadband heteronuclear decoupling of N15 during
		     acquisition
    dm(C13) = 'nn':  no broadband decoupling of C13 during acquisition
              'ny':  broadband heteronuclear decoupling of C13 during
                     acquisition

      phase = 1,2:  hypercomplex experiment with F1 quadrature (complex F1-FT)
     phase2 = 1,2:  hypercomplex experiment with F2 quadrature (complex F2-FT)
     phase3 = 1,2:  hypercomplex experiment with F3 quadrature (complex F3-FT)



   Recommended Values:  (MW = 16 kD)

	  jca = 140 Hz
	tauch = 3.57 ms
          jcn = 20-24 Hz
          jnh = 105 Hz
        taunh = 5.55 ms
*/

#include <standard.h>
#include <math.h>

#define MIN_DELAY	0.2e-6	/* shortest executable delay	*/
#define MIN_J		0.1	/* Hz				*/
#define MAX_SSTRIM	0.1	/* sec				*/


static int	phs1[1] = {0},
		phs2[2] = {0,2},
		phs3[2] = {1,3},
		phs4[4] = {0,1,2,3},
		rec[16] = {0,2,2,0,2,0,0,2,2,0,0,2,0,2,2,0};

static double	d2_init = 0.0,
		d3_init = 0.0,
		d4_init = 0.0;


/*-------------------------------
|				|
|       pulsesequence()/0	|
|				|
+------------------------------*/
pulsesequence()
{
/* Variable declaration */
   char		satflg[MAXSTR],
                sspul[MAXSTR],
                coshape[MAXSTR],
                fad1[MAXSTR],
                fad2[MAXSTR],
                fad3[MAXSTR],
                f1180[MAXSTR],
                f2180[MAXSTR],
                f3180[MAXSTR],
                prtinfo[MAXSTR],
                hdshape[MAXSTR],
                *dmc13,
                *dmn15;
   int          phase,
                phase2,
                phase3,
                cocapwrmove,
                satmove,
                t1_counter,
                t2_counter,
                t3_counter,
                c13dev,
                n15dev;
   double       ss,
                sstrim,
                sw1,
                sw2,
                sw3,
		t1evol,
		t2evol,
		t3evol,
		delayval,
                pwc13,
                pwc13lvl,
                pwca,
                pwcalvl,
		pwco,
		pwcolvl,
                cafrq,
                pwn15lvl,
                pwn15,
                satfrq,
                satdly,
                scuba,
                satpwr,
		hd90,
		hdlvl,
                hdres,
		jch,
		deltach,
		tauch,
		jcn,
		deltacn,
		jnh,
		deltanh,
		taunh;


/* Load variables */
   satfrq = getval("satfrq");
   satdly = getval("satdly");
   scuba = getval("scuba");
   satpwr = getval("satpwr");
   hdlvl = getval("hdlvl");
   hdres = getval("hdres");
   hd90 = getval("hd90");
   jch = getval("jch");
   deltach = getval("deltach");
   tauch = getval("tauch");
   jcn = getval("jcn");
   jnh = getval("jnh");
   deltanh = getval("deltanh");
   taunh = getval("taunh");
   pwca = getval("pwca");
   pwcalvl = getval("pwcalvl");
   pwco = getval("pwco");
   pwcolvl = getval("pwcolvl");
   pwc13 = getval("pwc13");
   pwc13lvl = getval("pwc13lvl");
   cafrq = getval("cafrq");
   pwn15lvl = getval("pwn15lvl");
   pwn15 = getval("pwn15");
   ss = getval("ss");
   sw1 = getval("sw1");
   sw2 = getval("sw2");
   sw3 = getval("sw3");
   sstrim = getval("sstrim");
   phase = (int) (getval("phase") + 0.5);
   phase2 = (int) (getval("phase2") + 0.5);
   phase3 = (int) (getval("phase3") + 0.5);
   c13dev = (int) (getval("c13dev") + 0.5);
   n15dev = (int) (getval("n15dev") + 0.5);

   getstr("coshape", coshape);
   getstr("hdshape", hdshape);
   getstr("sspul", sspul);
   getstr("fad1", fad1);
   getstr("fad2", fad2);
   getstr("fad3", fad3);
   getstr("satflg", satflg);
   getstr("f1180", f1180);
   getstr("f2180", f2180);
   getstr("f3180", f3180);
   getstr("prtinfo", prtinfo);


/* Initialize variables */
   satmove = ( fabs(tof - satfrq) >= 0.1 );
   cocapwrmove = ( fabs(pwcolvl - pwcalvl) >= 0.99 );

   deltacn = ( (jcn > MIN_J) ? 1/(2*jcn) : 0.0 );
   if (tauch < MIN_DELAY)
      tauch = ( (jch > MIN_J) ? 1/(2*jch) : 0.0 );
   if (deltach < MIN_DELAY)
      deltach = ( (jch > MIN_J) ? 1/(2*jch) : 0.0 );
   if (taunh < MIN_DELAY)
      taunh = ( (jnh > MIN_J) ? 1/(2*jnh) : 0.0 );
   if (deltanh < MIN_DELAY)
      deltanh = ( (jnh > MIN_J) ? 1/(2*jnh) : 0.0 );


/* Load phase cycles */
   settable(t1, 1, phs1);
   settable(t2, 2, phs2);
   settable(t4, 2, phs2);
   settable(t6, 2, phs2);
   settable(t3, 2, phs3);
   settable(t5, 4, phs4);
   settable(t8, 16, rec);

   setdivnfactor(t4, 2);
   setdivnfactor(t6, 4);
   setdivnfactor(t5, 8);


/* Check conditions */
   if ( (n15dev != DODEV) && (n15dev != DO2DEV) )
   {
      text_error("\ninvalid RF device for N15\n");
      abort(1);
   }
   else if ( (c13dev != DODEV) && (c13dev != DO2DEV) )
   {
      text_error("\ninvalid RF device for C13\n");
      abort(1);
   }
   else if (c13dev == n15dev)
   {
      text_error("\nN15 and C13 RF devices must be different\n");
      abort(1);
   }

   dmn15 = ( (n15dev == DODEV) ? dm : dm2 );
   dmc13 = ( (c13dev == DODEV) ? dm : dm2 );

   if ( (strcmp(dmn15, "n") != 0) && (strcmp(dmn15, "nn") != 0) &&
        (strcmp(dmn15, "ny") != 0) )
   {
      text_error("\n`dm(N15)` must be to 'n', 'nn', or 'ny'\n");
      abort(1);
   }

   if ( (strcmp(dmc13, "n") != 0) && (strcmp(dmc13, "nn") != 0) )
   {
      text_error("\n`dm(C13)` must be to 'n' or 'nn'\n");
      abort(1);
   }
 
   if (sstrim > MAX_SSTRIM)
   {
      text_error("\n`sstrim` is > maximum value\n");
      abort(1);
   }

   if ( fabs( cafrq - ((c13dev == DODEV) ? dof : dof2) ) >= 0.1 )
   {
      text_error("\ndof(C13) should be set to cafrq");
      abort(1);
   }

   if ( taunh > (deltacn/2) )
   {
      text_error("\ndeltacn > 2*taunh is a required condition");
      abort(1);
   }
   else if ( tauch > (deltacn/2) )
   {
      text_error("\ndeltacn > 2*tauch is a required condition");
      abort(1);
   }

   if (ix == 1)
   {
      if (prtinfo[A] == 'y')
      {
         double	maxt2,
		maxt3;

         maxt2 = deltacn - 2*tauch - 4.0e-6;
         t2_counter = (int) ( (maxt2 - d3_init)*sw2 + 0.5 );
         (void) printf("\nMaximum Ca t2 evolution time = %7.3lf ms\n",
			   maxt2*1e+3);
         (void) printf("Maximum Ca t2 increment = %4.0lf\n",
			   (double)t2_counter);

         maxt3 = deltacn - 2*POWER_DELAY;
         if (pwc13 > pwn15) 
            maxt3 -= (pwc13 - pwn15)/2;
         t3_counter = (int) ( (maxt3 - d4_init)*sw3 + 0.5 );
         (void) printf("Maximum N15 t3 evolution time = %7.3lf ms\n",
			   maxt3*1e+3);
         (void) printf("Maximum N15 t3 increment = %4.0lf\n\n",
			   (double)t3_counter);
      }
   }


/* Determine steady-state mode */
   if (ss < 0)
   {
      ss *= (-1);
      initval(ss, ssval);
      initval(ss, ssctr);
   }


/* Add in States-Haberkorn elements */
   if (phase == 2)		/* Ha t1 element	*/
      tsadd(t1, 1, 4);

   if (phase2 == 2)		/* Ca t2 element	*/
      tssub(t4, 1, 4);

   if (phase3 == 2)		/* N t3 element		*/
      tssub(t6, 1, 4);


/* Add in FAD */
   if (fad1[A] == 'y')
   {
      if (ix == 1)
        d2_init = d2;

      t1_counter = (int) ( (d2 - d2_init)*sw1 + 0.5 );
      if (t1_counter % 2)
      {
         tsadd(t1, 2, 4);       /* first Ha 90-degree pulse     */
         tsadd(t8, 2, 4);       /* receiver phase cycle         */
      }
   }

   if (fad2[A] == 'y')
   {
      if (ix == 1) 
        d3_init = d3; 
 
      t2_counter = (int) ( (d3 - d3_init)*sw2 + 0.5 ); 
      if (t2_counter % 2) 
      { 
         tsadd(t4, 2, 4);       /* first Ca 90-degree pulse     */ 
         tsadd(t8, 2, 4);       /* receiver phase cycle         */ 
      } 
   }

   if (fad3[A] == 'y')
   {
      if (ix == 1) 
        d4_init = d4; 
 
      t3_counter = (int) ( (d4 - d4_init)*sw3 + 0.5 ); 
      if (t3_counter % 2) 
      { 
         tsadd(t6, 2, 4);       /* first N15 90-degree pulse    */ 
         tsadd(t8, 2, 4);       /* receiver phase cycle         */ 
      } 
   }


/* BEGIN ACTUAL PULSE SEQUENCE CODE */
   status(A);
      rlpower(tpwr, TODEV);		/* H1 hard-pulse power level	*/
      rlpower(pwn15lvl, n15dev);	/* N15 hard-pulse power level	*/
      rlpower(pwc13lvl, c13dev);	/* C13 hard-pulse power level	*/

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
      if (satflg[A] == 'y')
      {
         if (satmove)
            offset(satfrq, TODEV);

         rlpower(satpwr, TODEV);
         rgpulse(satdly, zero, 4.0e-5, 0.2e-6);
         if (satmove)
            offset(tof, TODEV);

         rlpower(tpwr, TODEV);
         rcvroff();
         delay(4.0e-5);

	 if (scuba > MIN_DELAY)
         {
	    rgpulse(pw, zero, 1.0e-6, 0.0);
	    rgpulse(2.0*pw, one, 1.0e-6, 0.0);
	    rgpulse(pw, zero, 1.0e-6, 0.0);
	    delay(scuba);
	    rgpulse(pw, zero, 1.0e-6, 0.0);
	    rgpulse(2.0*pw, one, 1.0e-6, 0.0);
	    rgpulse(pw, zero, 1.0e-6, 0.0);
	    delay(scuba);
	 }
      }


/* Start the pulse train */
      rcvroff();
      rgpulse(pw, t1, rof1, 0.0);


/* Ha evolution time and INEPT transfer from Ha to Ca */
      t1evol = d2/2;
      if (t1evol > MIN_DELAY)
      {
         txphase(zero);
         delay( t1evol + (deltach/2) - 2*pwc13 - 2.0e-6 );
         genrgpulse(2*pwc13, t2, 2.0e-6, 0.0, c13dev);
         delay(t1evol);
         rgpulse(2*pw, zero, 0.0, 0.0);
      }
      else
      {
         delay( (deltach/2) - 2.0e-6 );
         gensim2pulse(2*pw, 2*pwc13, zero, t2, 2.0e-6, 0.0, TODEV, c13dev);
      }

      delay( (deltach/2) - 2.0e-6 );
      gensim2pulse(pw, pwc13, t3, t4, 2.0e-6, 0.0, TODEV, c13dev);


/* Refocus Ha to Ca INEPT transfer and start INEPT transfer from Ca to N15;
   constant-time Ca evolution time */
      rlpower(pwcalvl + 6.0, c13dev);	/* Ca power level for 180	*/
      rlpower(hdlvl, TODEV);		/* H1 decoupling power level	*/
      obsprgon(hdshape, hd90, hdres);	/* start H1-channel WFG		*/
      txphase(zero);

      delayval = tauch - 2*POWER_DELAY - PRG_START_DELAY;
      if (pw > pwc13)
         delayval -= (pw - pwc13)/2;

      delay(delayval);
      xmtron();				/* turn on H1 RF for decoupling	*/


/* Continue INEPT transfer from Ca(t2) to N15 with CO decoupling */
      t2evol = d3/2;
      delayval = (deltacn/2) - tauch - t2evol - 2.0e-6;
      if (delayval < MIN_DELAY)
      {
         text_error("\nCa t2 evolution time is too large for CT");
         abort(1);
      }

      delay(delayval);
      gensim2pulse(pwca, 2*pwn15, t5, zero, 2.0e-6, 0.0, c13dev, n15dev);

      delayval = t2evol - WFG_STOP_DELAY - POWER_DELAY - 1.0e-5;
      if (pwn15 > pwc13)
         delayval -= (pwn15 - pwc13)/2;

      if (delayval > MIN_DELAY)
      {
         if (cocapwrmove)
         {
            delay( (deltacn/2) - pwco - WFG_START_DELAY - POWER_DELAY
			- 12.0e-6 );
            rlpower(pwcolvl + 6.0, c13dev);
            delay(1.0e-5);
         }
         else
         {
            delay( (deltacn/2) - pwco - WFG_START_DELAY - 2.0e-6 );
         }

         genshaped_pulse(coshape, pwco, zero, 2.0e-6, 0.0, 0.0, 0.0, c13dev);

         delay(1.0e-5);
         rlpower(pwc13lvl, c13dev);	/* hard pulse C13 power level	*/
         delay(delayval);
      }
      else
      {
         delayval = (deltacn/2) + t2evol - pwco - POWER_DELAY -
			WFG_START_DELAY - WFG_STOP_DELAY - 1.2e-5;
         if (pwn15 > pwc13)
            delayval -= (pwn15 - pwc13)/2;

         if (cocapwrmove)
         {
            rlpower(pwcolvl + 6.0, c13dev);
            delayval -= POWER_DELAY;
         }

         delay(delayval);
         genshaped_pulse(coshape, pwco, zero, 2.0e-6, 0.0, 0.0, 0.0, c13dev);

         rlpower(pwc13lvl, c13dev);	/* hard pulse C13 power level	*/
         delay(1.0e-5);
      }

      gensim2pulse(pwc13, pwn15, zero, t6, 2.0e-6, 0.0, c13dev, n15dev);


/* Refocus Ca-->N15 INEPT transfer during CT N15(t3) with CO decoupling */
      delay(1.0e-5);
      rlpower(pwcalvl + 6.0, c13dev);	/* Ca power level for sel. 180	*/

      t3evol = d4/2;
      delayval = (deltacn/2) - t3evol - POWER_DELAY - 1.0e-5;
      if (pwc13 > pwn15)
         delayval -= (pwc13 - pwn15)/2;

      if (delayval < MIN_DELAY)
      {
         text_error("\nN15 t3 evolution time is too large for CT");
         abort(1);
      }

      delay(delayval);
      gensim2pulse(pwca, 2*pwn15, zero, zero, 0.0, 0.0, c13dev, n15dev);

      if ( t3evol > (taunh + WFG_STOP_DELAY + MIN_DELAY) )
      {
         if (cocapwrmove)
         {
            delay( (deltacn/2) - pwco - WFG_START_DELAY - POWER_DELAY
			- 1.0e-5 );
            rlpower(pwcolvl + 6.0, c13dev);
            delay(1.0e-5);
         }
         else
         {
            delay( (deltacn/2) - pwco - WFG_START_DELAY );
         }

         genshaped_pulse(coshape, pwco, zero, 0.0, 0.0, 0.0, 0.0, c13dev);
         delay(t3evol - taunh - WFG_STOP_DELAY);
      }
      else
      {
         delay( (deltacn/2) - taunh + t3evol );
      }


/* Stop H1 decoupling */
      xmtroff();			/* turn off H1 RF for BBD       */
      rlpower(tpwr, TODEV);		/* H1 hard pulse power          */
      txphase(one);
      obsprgoff();			/* turn off H1 WFG for BBD      */
      rgpulse(pw, one, 0.0, 0.0);


/* Start INEPT transfer from N15 to Hn and finish CT t3 time */
      delayval = taunh - pw - POWER_DELAY - PRG_STOP_DELAY - 2.0e-6;
      if (pw > pwn15)
         delayval -= (pw - pwn15)/2;

      if ( t3evol < (taunh + WFG_STOP_DELAY + MIN_DELAY) )
      { /* Do CO 180 during taunh under certain conditions */
         delayval -= pwco + WFG_START_DELAY;
         if (cocapwrmove)
         {
            delayval -= POWER_DELAY + 1.0e-5;
            rlpower(pwcolvl + 6.0, c13dev);
            delay(1.0e-5);
         }

         if ( (delayval - t3evol) < MIN_DELAY )
         {
            genshaped_pulse(coshape, pwco, zero, 0.0, 0.0, 0.0, 0.0, c13dev);
            delay(delayval);
         }
         else
         {
            delay(delayval - t3evol);
            genshaped_pulse(coshape, pwco, zero, 0.0, 0.0, 0.0, 0.0, c13dev);
            delay(t3evol);
         }
      }
      else
      {
         delay(delayval);
      }

      gensim2pulse(pw, pwn15, zero, zero, 2.0e-6, 0.0, TODEV, n15dev);


/* Refocus INEPT transfer from N15 to Hn */
      delay(deltanh/2);
      gensim2pulse(2*pw, 2*pwn15, zero, zero, 0.0, 0.0, TODEV, n15dev);

      rlpower(dpwr, DODEV);		/* N15 or C13 decoupling power  */
      rlpower(dpwr2, DO2DEV);		/* C13 or N15 decoupling power  */
      genqdphase(zero, n15dev);		/* N15 decoupling phase		*/
      genqdphase(zero, c13dev);		/* C13 final phase		*/

      delay(rof2);
      rcvron();

      delay( (deltanh/2) - rof2 - 2*POWER_DELAY - PRG_START_DELAY );

    
/* Start N15 broadband decoupling */
   status(B);
      setreceiver(t8);
}
