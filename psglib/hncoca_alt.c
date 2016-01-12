/* hncoca_alt.c  -  modified HN(CO)CA 3D experiment by Ikura and Bax; uses
		    shifted lminar pulses for off-resonance excitation of
		    the CO spins; uses a broadband inversion pulse to refocus
		    the Ca-N15 and CO-N15 scalar interactions during t2; uses
		    broadband H1 decoupling during t1, t2, both tau(cn) periods,
		    and both tau(cc) periods; does not refocus the Ca-N15
		    scalar interactions during the two tau(cc) periods; the
		    C13 and N15 RF channels are user-selectable.

                        has two less 15N pulses than in hncoca.c . This may give
                     somewhat better sensitivity at the expense of less selectivity
                     (intra-residue connectivities may be seen)

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
      hdlvl = power level for H1 broadband decoupling during t(cn), t1, and t2
    hdshape = H1 broadband decoupling sequence
      hdres = tip-angle resolution for H1 broadband decoupling sequence
       hd90 = 90 degree pulse length for H1 broadband decoupling
      hdfrq = NOT DEFINED ==> tof is assumed to define the center of the H1
              broadband decoupling profile
         hs = 'yn':  homospoil pulse (hst) during the d1 relaxation delay
        jnh = one-bond H-N heteronuclear coupling constant (in Hz)
    deltanh = if deltanh < 0.2e-6, deltanh = 1/(2*jnh)
      taunh = should be set to 5.55 ms (or 1/[2*90]); if taunh < 0.2e-6,
              taunh = 1/(2*jnh)
        jcn = one-bond heteronuclear N-CO coupling constant (in Hz)
        jcc = one-bond CO-Ca coupling constant (in Hz)
       tpwr = power level for 1H transmitter pulses
         pw = 90 degree xmtr pulse length for protons (the observed nucleus)
      cafrq = frequency for Ca spins
      cofrq = frequency for CO spins
   dof(C13) = (cafrq + cofrq)/2 for maximum sensitivity
    pwcalvl = power level for Ca decoupler pulses
       pwca = 90 degree decoupler pulse length for Ca at `pwcalvl`
    coshape = frequency-shifted pulse for CO excitation
    n15shape= shape for N15 180
    pwcolvl = power level for 90-degree C0 decoupler pulses
       pwco = 90 degree decoupler pulse length for C0 at `pwcolvl`
   pwcoalvl = power level for hard C13 pulse
      pwcoa = 90 degree decoupler pulse length for C13 at `pwcoalvl`
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
                 modified from hncoca4FSP.c(S.Farmer,Varian)  1/27/92
*/



#include <standard.h>
#include <math.h>

#define MIN_DELAY	0.2e-6		/* shortest executable delay	*/
#define MIN_J		0.1		/* in Hz			*/
#define MIN_NULL        0.0001          /* sec                          */

#define WFG2_START_DELAY	20.2e-6
#define WFG2_STOP_DELAY		8.2e-6


static int	  count = 0,
		 phs[2] = {0,2},
		rec[64] = {0,2,2,0,0,2,2,0,2,0,0,2,2,0,0,2,
			   0,2,2,0,0,2,2,0,2,0,0,2,2,0,0,2,
			   2,0,0,2,2,0,0,2,0,2,2,0,0,2,2,0,
			   2,0,0,2,2,0,0,2,0,2,2,0,0,2,2,0};

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
		fad2[MAXSTR],
		f1180[MAXSTR],
		f2180[MAXSTR],
		hdshape[MAXSTR],
		coshape[MAXSTR],
                n15shape[MAXSTR],
		prtinfo[MAXSTR],
		*dmc13,
		*dmn15;
   int          phase,
		phase2,
		satmove,
		exp3D,
		t1_counter,
		t2_counter,
		c13dev,
		n15dev;
   double       ss,
 		delayval,
		d2local,
		t1evol_1,
		t1evol_2,
		sw1,
		sw2,
		cafrq,
                pwcalvl,
                pwca,
                pwcolvl,
                pwco,
                pwcoalvl,
                pwcoa,
                pwn15lvl,
                pwn15,
                satfrq,
                satdly,
                satpwr,
		hdlvl,
		hdres,
		hd90,
		jnh,
		deltanh,
		taunh,
		jcn,
		deltacn,
		jcc,
		deltacc;


/* Load variables */
   satfrq = getval("satfrq");
   satdly = getval("satdly");
   satpwr = getval("satpwr");
   cafrq = getval("cafrq");
   pwcalvl = getval("pwcalvl");
   pwca = getval("pwca");
   pwcolvl = getval("pwcolvl");
   pwco = getval("pwco");
   pwcoalvl = getval("pwcoalvl");
   pwcoa = getval("pwcoa");
   pwn15lvl = getval("pwn15lvl");
   pwn15 = getval("pwn15");
   hd90 = getval("hd90");
   hdlvl = getval("hdlvl");
   hdres = getval("hdres");
   taunh = getval("taunh");
   deltanh = getval("deltanh");
   jnh = getval("jnh");
   jcn = getval("jcn");
   jcc = getval("jcc");
   ss = getval("ss");
   sw1 = getval("sw1");
   sw2 = getval("sw2");
   phase = (int) (getval("phase") + 0.5);
   phase2 = (int) (getval("phase2") + 0.5);
   c13dev = (int) (getval("c13dev") + 0.5);
   n15dev = (int) (getval("n15dev") + 0.5);

   getstr("hdshape", hdshape);
   getstr("coshape", coshape);
   getstr("n15shape", n15shape);
   getstr("sspul", sspul);
   getstr("fad1", fad1);
   getstr("fad2", fad2);
   getstr("satmode", satmode);
   getstr("f1180", f1180);
   getstr("f2180", f2180);
   getstr("prtinfo", prtinfo);


/* INITIALIZE VARIABLES */
   deltacc = ( (jcc > MIN_J) ? 1/(2*jcc) : 0.0 );
   deltacn = ( (jcn > MIN_J) ? 1/(2*jcn) : 0.0 );
   if (taunh < MIN_DELAY)
      taunh = ( (jnh > MIN_J) ? 1/(2*jnh) : 0.0 );
   if (deltanh < MIN_DELAY)
      deltanh = ( (jnh > MIN_J) ? 1/(2*jnh) : 0.0 );



/* Load phase cycles */
   settable(t1, 2, phs);	/* Ca 90	*/
   settable(t2, 2, phs);	/* Ca 90	*/
   settable(t3, 2, phs);	/* CO 180	*/
   settable(t4, 2, phs);	/* H 90		*/
   settable(t5, 2, phs);	/* H 180	*/
   settable(t7, 2, phs);	/* N 180	*/
   settable(t9, 2, phs);	/* N 90		*/
   settable(t10, 64, rec);

   setdivnfactor(t9, 2);
   setdivnfactor(t3, 4);
   setdivnfactor(t2, 8);
   setdivnfactor(t5, 16);
   setdivnfactor(t4, 32);
   setdivnfactor(t7, 64);


/* Check conditions */
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


/* Determine steady-state mode */
   if (ss < 0)
   {
      ss *= (-1);
      initval(ss, ssval);
      initval(ss, ssctr);
   }


/* Add in States-Haberkorn element */
   if (phase == 2)		/* C13 t1 element */
      tsadd(t1, 1, 4);

   assign(zero, v1);
   if (phase2 == 2)		/* N15 t2 element */
      sub(v1, one, v1);


/* Add in FAD */
   if (fad1[A] == 'y')		/* for C13 */
   {
      if (ix == 1)
        d2_init = d2;
 
      t1_counter = (int) ( (d2 - d2_init)*sw1 + 0.5 );
      if (t1_counter % 2)
      {   
         tsadd(t1, 2, 4);       /* first C13a 90-degree pulse   */
         tsadd(t10, 2, 4);      /* receiver phase cycle         */
      }
   }

   if (fad2[A] == 'y')		/* for N15 */
   {
      if (ix == 1)
        d3_init = d3;

      t2_counter = (int) ( (d3 - d3_init)*sw2 + 0.5 );
      if (t2_counter % 2)
      {
         add(two, v1, v1);	/* first N15 90-degree pulse    */
         tsadd(t10, 2, 4);	/* receiver phase cycle		*/
      }
   }


/* BEGIN ACTUAL PULSE SEQUENCE CODE */
   status(A);
      rlpower(tpwr, TODEV);		/* H1 hard-pulse power level	*/
      rlpower(pwn15lvl, n15dev);	/* N15 hard-pulse power level	*/
      rlpower(pwcalvl + 6.0, c13dev); 	/* Ca selective 180 power level	*/
      offset(cafrq, c13dev);            /* Ca freq                      */

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
      genqdphase(zero, c13dev);
      rgpulse(pw, zero, rof1, 5.0e-6);

      txphase(zero);
      genqdphase(v1, n15dev);
      delay(deltanh - 5.0e-6);
      gensim2pulse(2*pw, pwn15, zero, v1, 0.0, 5.0e-6, TODEV, n15dev);

      txphase(one);
      delay(deltanh - 5.0e-6);
      rgpulse(pw, one, 0.0, 0.0);

      delay( (deltacn/2) - deltanh );
      genrgpulse(pwca, zero, 0.0, 0.0, c13dev);
		/* 180-degree pulse due to 6 dB more power */

      delayval = taunh + deltanh - (deltacn/2);
      delay(delayval - 2*POWER_DELAY - PRG_START_DELAY);


/* start H1 decoupling now! */
      txphase(zero);			/* H1 phase for BB decouple	*/
      rlpower(hdlvl, TODEV);		/* power for H1 BB decouple	*/
      rlpower(pwcolvl, c13dev);
      obsprgon(hdshape, hd90, hdres);   /* start up H1 WFG for BBD      */
      xmtron();				/* turn on H1 RF		*/


/* Establish HMQC between N15 and CO */
      genqdphase(zero, c13dev);
      delay(deltacn - taunh - deltanh - WFG_START_DELAY);
      genshaped_pulse(coshape, pwco, zero, 0.0, 0.0, 0.0, 0.0, c13dev);


/* Establish HMQC between CO and Ca */
      rlpower(pwcalvl, c13dev);
      genqdphase(t1, c13dev);
      delay( (deltacc/2) - POWER_DELAY - WFG_STOP_DELAY );
      genqdphase(t7, n15dev);
      delay(deltacc/2);


/* Adjust t1 times for proper F1 phasing */
      if (ix == 1)
      {
         t1adj_1 = (2*pwca)/M_PI + (pwco/4) + POWER_DELAY + WFG2_START_DELAY;
         t1adj_2 = (2*pwca)/M_PI + (pwco/4) + POWER_DELAY + WFG2_STOP_DELAY;
         t1evol_1 = (4 * sw1) * ( MIN_DELAY + t1adj_1 );
         t1evol_2 = (4 * sw1) * ( MIN_DELAY + t1adj_2 );

         tmpd2 = ( (t1evol_1 > t1evol_2) ? t1evol_1 : t1evol_2 );
         count = (int) (tmpd2) + 1;

         if (prtinfo[0] == 'y')
         {
            char   msge[128];
            double maxsw1;

            maxsw1 = ( (f1180[0] == 'y') ? (sw1/tmpd2) : (2*sw1)/tmpd2 );
            (void) sprintf(msge, "\nMaximum C13 spectral width = %d\n",
                                (int) (maxsw1 + 0.5) );
            text_error(msge);
         }
 
         tmpd2 = count * (0.5/sw1);
      }

      if ( (f1180[0] != 'y') && (count < 3) )
      {
         t1evol_1 = (d2/2) - t1adj_1;
         t1evol_2 = (d2/2) - t1adj_2;
      }
      else
      {   
         if ( (count > 2) || ((f1180[0] == 'y') && (count > 1)) )
         {
            if (ix == 1)
               text_error("WARNING:  sw1 is too large for proper t1 timing\n");
         }

         d2local = d2 + tmpd2;
         t1evol_1 = (d2local/2) - t1adj_1;
         t1evol_2 = (d2local/2) - t1adj_2;
      }

      if (t1evol_1 < MIN_DELAY)
         t1evol_1 = 0.0;
      if (t1evol_2 < MIN_DELAY)
         t1evol_2 = 0.0;


/* Start Ca evolution time after 90-degree pulse */
      genrgpulse(pwca, t1, 0.0, 0.0, c13dev);
      genqdphase(t3, c13dev);
      rlpower(pwcolvl + 6.0, c13dev);	/* CO selective 180 power level	*/
      delay(t1evol_1);
      gensim2shaped_pulse(n15shape, coshape, 2*pwn15, pwco, t7, t3, 0.0,
				0.0, 0.0, 0.0, n15dev, c13dev);
		/* 180-degree pulse on CO due to 6 dB more power */

      rlpower(pwcalvl, c13dev);
      genqdphase(t2, c13dev);
      delay(t1evol_2);

      genrgpulse(pwca, t2, 0.0, 0.0, c13dev);
      genqdphase(zero, c13dev);


/* Remove HMQC between Ca and CO */
      delay(deltacc/2);
      genqdphase(t9, n15dev);
      rlpower(pwcolvl, c13dev);

      delay( (deltacc/2) - POWER_DELAY - WFG_START_DELAY );
      genshaped_pulse(coshape, pwco, zero, 0.0, 0.0, 0.0, 0.0, c13dev);


/* Start to remove HMQC between CO and N15 */
      rlpower(pwcoalvl, c13dev);
      genqdphase(zero, c13dev);
      delayval = (deltacn - taunh - deltanh)/2 - POWER_DELAY - WFG_STOP_DELAY;

      if (exp3D)
      {
         delayval -= (8*pwcoa) + ((c13dev == DO2DEV) ? OFFSET_LTCH_DELAY :
			OFFSET_DELAY);
         offset( ((c13dev == DO2DEV) ? dof2 : dof), c13dev);
      }

      delay(delayval);


/* 3D section */
      if (exp3D)
      {
         delay(d3/2);	  /* N15 evolution time:  function of ni2 and sw2 */
         genrgpulse(pwcoa, zero, 0.0, 0.0, c13dev);
         genrgpulse(2*pwcoa, one, 0.0, 0.0, c13dev);
         genrgpulse(pwcoa, zero, 0.0, 0.0, c13dev);
         delay(d3/2);
         genrgpulse(pwcoa, zero, 0.0, 0.0, c13dev);
         genrgpulse(2*pwcoa, one, 0.0, 0.0, c13dev);
         genrgpulse(pwcoa, zero, 0.0, 0.0, c13dev);
      }


/* Remove more of the HMQC between CO and N15 */
      delayval = (deltacn - taunh - deltanh)/2 - pw - POWER_DELAY -
			PRG_STOP_DELAY;
      delay(delayval);


/* Stop H1 decoupling */
      xmtroff();                        /* turn off H1 RF for BBD       */
      rlpower(tpwr, TODEV);
      txphase(one);
      obsprgoff();                      /* turn of H1 WFG for BBD       */
      rgpulse(pw, one, 0.0, 0.0);	/* H1 purge pulse               */


/* Establish NxHz antiphase state and complete refocusing of 
   NxCOz antiphase state */
      genqdphase(zero, c13dev);
      delayval = taunh + deltanh - (deltacn/2) - POWER_DELAY;

      if (exp3D)
      {
         delayval -= ((c13dev == DO2DEV) ? OFFSET_LTCH_DELAY :
				OFFSET_DELAY);
         offset(cafrq, c13dev);
      }

      rlpower(pwcalvl + 6.0, c13dev);	/* Ca selective 180 pulse power	*/
      delay(delayval);
      genrgpulse(pwca, zero, 0.0, 0.0, c13dev);
			/* 180-degree pulse due to 6 dB more power	*/

      txphase(t4);
      delay( (deltacn/2) - deltanh );
      rgpulse(pw, t4, 0.0, 5.0e-6);
      txphase(t5);
      delay(deltanh - 5.0e-6);
      gensim2pulse(2*pw, pwn15, t5, t9, 0.0, 0.0, TODEV, n15dev);

      rlpower(dpwr, DODEV);             /* N15 or C13 decoupling power */
      rlpower(dpwr2, DO2DEV);           /* C13 or N15 decoupling power */
      genqdphase(zero, n15dev);         /* N15 decoupling phase */
      genqdphase(zero, c13dev);         /* C13 final phase */
      delay(rof2);
      rcvron();

      delay(deltanh - 2*POWER_DELAY - PRG_START_DELAY - rof2);


/* Start N15 broadband decoupling */
   status(B);
      setreceiver(t10);
}
