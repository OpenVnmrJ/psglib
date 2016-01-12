/* hcacon_ct.c      -  HCA(CO)N 3D experiment by Fairbrother et al.; uses a
	   	    frequency-shifted pulse on the CO spins
                    with modifications to allow 1H 180 to cross CO 90 at long t1 times
                    sw1, ni and phase for Ca (t1)
                    sw2, ni2 and phase2 for 15N (t2)
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
   pwc13lvl = first power level for Ca decoupler pulses
      pwc13 = 90 degree decoupler pulse length for Ca at `pwc13lvl`
   pwn15lvl = power level for N15 decoupler pulses
      pwn15 = 90 degree decoupler pulse length for N15
       dpwr = power level for C13 broadband decoupling
        dof = decoupler frequency-should be centered on alpha carbons
        jch = CH scalar coupling - used to calculate tau
       Del1 = second delay in sequence
       Del2 = third delay in sequence
        del = little delta delay
          ConstT =  constant time delay
    coshape = SLP pulse of 0 phase for off resonance CO
   cshape1 = SLP pulse of Ca phase 1 CO phase 0
   cshape2 = SLP pulse of Ca phase 1 CO phase 2
   cshape3 = SLP pulse Ca phase 0 CO phase 0, but lower power for Ca
         dm = 'nn':  no broadband decoupling of C13 during acquisition
              'ny':  broadband heteronuclear decoupling of C13
		        during acquisition
        dm2 = 'nn':  no broadband decoupling of N15 during acquisition

      phase = 1,2:  hypercomplex experiment with F1 quadrature (complex F1-FT)
     phase2 = 1,2:  hypercomplex experiment with F2 quadrature (complex F2-FT)
*/



#include <standard.h>
#include <math.h>

#define MIN_DELAY	0.2e-6		/* shortest executable delay */



static int	 phs1[1] = {0},
		 phs2[8] = {1,1,1,1,3,3,3,3},
		 phs3[2] = {0,2},
                 phs4[16]= {0,0,0,0,0,0,0,0,2,2,2,2,2,2,2,2},
                 phs5[4] = {0,0,2,2},
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
		coshape[MAXSTR],
		cshape1[MAXSTR],
		cshape2[MAXSTR],
                cshape3[MAXSTR];

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
                pwc13lvl,
                pwc13,
                pwn15lvl,
                pwn15,
                satfrq,
                satdly,
                satpwr,
		jch,
                tau,tau_,
                Del1,
                Del2,
                del,
                ConstT,
                del_,
                ttime,
                Del1_,
                Del2_,
                t1evol_,
                ConstT_;


/* Load variables */
   satfrq = getval("satfrq");
   satdly = getval("satdly");
   satpwr = getval("satpwr");
   pwc13lvl = getval("pwc13lvl");
   pwc13 = getval("pwc13");
   pwn15lvl = getval("pwn15lvl");
   pwn15 = getval("pwn15");
   jch = getval("jch");
   Del1 = getval("Del1");
   Del2 = getval("Del2");
   ConstT = getval("ConstT");
   del = ConstT - Del1;
   ss = getval("ss");
   sw1 = getval("sw1");
   sw2 = getval("sw2");
   phase = (int) (getval("phase") + 0.5);
   phase2 = (int) (getval("phase2") + 0.5);

   getstr("coshape", coshape);
   getstr("cshape1", cshape1);
   getstr("cshape2", cshape2);
   getstr("cshape3", cshape3);
   getstr("sspul", sspul);
   getstr("fad1", fad1);
   getstr("fad2", fad2);
   getstr("satmode", satmode);
   getstr("f1180", f1180);
   getstr("f2180", f2180);


/* Load phase cycles */
   settable(t1, 1, phs1); /* phase of 1st Ca 90 */
   settable(t2, 8, phs2); /* phase of 2nd Ca 90 */
   settable(t3, 2, phs3);/* phase of 1st N15 90 */ 
   settable(t4, 16, phs4);/* phase of 2nd N15 90 */ 
   settable(t5, 4, phs5); /* phase of 2nd CO 90 */
   settable(t6, 16, rec); /* phase of receiver */

   tau = 1/(4*jch); 

/* Check conditions */
   satmove = ( fabs(tof - satfrq) >= 0.1 );

   if ( (dm2[A] == 'y') || (dm2[B] == 'y') || (dm2[C] == 'y'))
   {
      text_error("`dm2` must be  'nn'  for N15\n");
      abort(1);
   }

   if ( (dm[A] == 'y') || (dm[B] == 'y') )
   {
      text_error("`dm` must be 'nnn' or 'nny' for C13\n");
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
   if (phase == 2)		/* Ca t1 element */
      tsadd(t1,1,4);

   if (phase2 == 2)		/* N15 t2 element */
      tsadd(t3, 1, 4);


/* Add in FAD */
   if (fad1[A] == 'y')		/* for Ca */
   {
      if (ix == 1)
        d2_init = d2;

      t1_counter = (int) ( (d2 - d2_init)*sw1 + 0.5 );
      if (t1_counter % 2)
      {
         tsadd(t1, 2, 4);	/* first Ca 90-degree pulse 	*/
         tsadd(t6, 2, 4);	/* receiver phase cycle		*/
      }
   }

   if (fad2[A] == 'y')		/* for 15N */
   {
      if (ix == 1)
        d3_init = d3;
 
      t2_counter = (int) ( (d3 - d3_init)*sw2 + 0.5 );
      if (t2_counter % 2)
      {   
         tsadd(t3, 2, 4);       /* first N15 90-degree pulse   */
         tsadd(t6, 2, 4);       /* receiver phase cycle         */
      }
   }


/* BEGIN ACTUAL PULSE SEQUENCE CODE */
   status(A);
      rlpower(tpwr, TODEV);		/*  H1 hard-pulse power level	*/
      rlpower(pwn15lvl, n15dev);	/* N15 hard-pulse power level	*/
      rlpower(pwc13lvl, c13dev); 	/* C13 pulse power		*/

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
      }

      rcvroff();

status(B);  
/* INEPT transfer from H to C */
      
      rgpulse(pw, zero, rof1, 0.0);
      
      decphase(zero);		/* C13 phase */
      delay(tau);

      simpulse(2*pw, 2*pwc13, zero, zero, 0.0, 0.0);
      
      txphase(one);			/* H1 phase  */
      decphase(t1);		/* C13 phase */
      delay(tau); 
      simpulse(pw, pwc13, one, t1, 0.0, 0.0);

/* Ca evol. is governed by ni and sw */
      t1evol = d2/2;
      if (f1180[0] == 'y')
         t1evol += 0.25/sw;

/* adjust the Del1 delay so that it will not = t1evol  */
      Del1 = (double) ((int) (Del1*sw1*2 + 0.5))/(sw1*2);
      Del1 +=   2e-6 + 2*pw; 


      txphase(zero);
      ttime=0.0;
      if (t1evol < Del1)       /* if statement allows t1/2 to exceed Del1 */
      {
        delay(t1evol);
        rgpulse(pw*2, zero, 0.0, 0.0); /* H1 180 pulse */
        ttime = t1evol + 2*pw;
        Del1_ = Del1 - ttime;
        delay(Del1_);
        decshaped_pulse(coshape,pwc13,zero,0.0,0.0);
        ttime += Del1_ + pwc13 + WFG_START_DELAY + WFG_STOP_DELAY;
      }
     else
      {
        delay(Del1);
        decshaped_pulse(coshape,pwc13,zero,0.0,0.0);
        ttime = Del1 + pwc13 + WFG_START_DELAY + WFG_STOP_DELAY;
        t1evol_ = t1evol - ttime;
        delay(t1evol_);
        rgpulse(pw*2, zero, 0.0, 0.0); /* H1 180 pulse */
        ttime += t1evol_ + 2*pw;
      }

      rlpower(pwc13lvl+6,c13dev);
      ConstT_ = ConstT/2 + t1evol - ttime - pwc13/2 - POWER_DELAY;
      decphase(zero);
      delay(ConstT_);
      decrgpulse(pwc13,zero,0.0,0.0);
      rlpower(pwc13lvl,c13dev);
      decphase(t2);
      dec2phase(t3);
      ConstT_ = ConstT/2 - t1evol - pwc13/2 - POWER_DELAY;
      delay(ConstT_);

/* 15N pulses and 15N evolution time with simultaneous Ca pulse */      

      sim3pulse(0.0,pwc13,pwn15,zero,t2,t3,0.0,0.0); 

      t2evol_1 = t2evol_2 = d3/2;	/* N15 is governed by ni2 and sw2 */
      if (f2180[0] == 'y')
      {
         t2evol_1 += 0.25/sw2;
         t2evol_2 += 0.25/sw2;
      }

      if ( (t2evol_1 > MIN_DELAY) || (t2evol_2 > MIN_DELAY) )
      {
         t2evol_1 -= (2*pwn15/M_PI) + pwc13/2 + POWER_DELAY + WFG_START_DELAY;
         t2evol_2 -= (2*pwn15/M_PI) + pwc13/2 + WFG_STOP_DELAY;

         if (t2evol_1 < MIN_DELAY)
            t2evol_1 = 0.0;
         if (t2evol_2 < MIN_DELAY)
            t2evol_2 = 0.0;
       }

      
      decphase(zero);
      rlpower(pwc13lvl+6,c13dev);
     
      delay(t2evol_1);

/* 180 on CO */
      decshaped_pulse(coshape, pwc13, zero, 0.0, 0.0); 

      dec2phase(t4);
      delay(t2evol_2);
      dec2rgpulse(pwn15, t4, 0.0, 0.0); /* 2nd N15 90 pulse */

      decphase(zero);
      delay(del/2);

      decrgpulse(pwc13,zero,0.0,0.0); /*  180 on Ca */

      decphase(zero);
      del_ = del/2 - WFG_START_DELAY;
      delay(del_);

/* Simultaneous 90 pulses on Ca and CO */
/* C power already 6 db higher to account for power split */

      getelem(t5,ct,v2);
      ifzero(v2);
      decshaped_pulse(cshape1,pwc13,zero,0.0,0.0);
      elsenz(v2);
      decshaped_pulse(cshape2,pwc13,zero,0.0,0.0);
      endif(v2);

      Del2_ = Del2/2 - POWER_DELAY - WFG_STOP_DELAY - WFG_START_DELAY;

      delay(Del2_);

      decshaped_pulse(cshape3,4*pwc13, zero, 0.0, 0.0);  /* simultaneous C 180s */
      rlpower(pwc13lvl, DODEV);
      Del2_ = Del2/2  - WFG_STOP_DELAY;
      delay(Del2_);

      simpulse(pw, pwc13, zero, zero, 0.0, 0.0);

      delay(tau);

      simpulse(2*pw, 2*pwc13, zero, zero, 0.0, rof2);
      rlpower(dpwr, DODEV);
      rlpower(dpwr2, DO2DEV);
      tau_=tau-2*POWER_DELAY;
      delay(tau_);

      rcvron();
   status(C);
      setreceiver(t6);
}
