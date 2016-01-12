/* hmqcjr_3rf- heteronuclear multiple-quantum coherence with hypercomplex
         implementation, steady-state transients for either all or
         just the first t1 increment, pulse+receiver phasecycling during
         steady-state transients, proton pulses are jump-return pulses

          USES THE THIRD CHANNEL FOR X PULSES AND DECOUPLING

   ref: Summers et al., j. amer. chem. soc. 108:4285-4294 (1986)


   Parameters:

        dof2 = frequency for decoupling during acquisition
         hs = 'yn':  homospoil pulse (hst) during the d1 relaxation delay
         pw = 90 degree xmtr pulse length for protons (the observed nucleus)
      p1 = correction to second pulse of jr pair 
     pwx2lvl = power level for decoupler pulses
        pwx2 = 90 degree decoupler pulse length for the heteronucleus
          j = one-bond heteronuclear coupling constant (in Hz)
      tau =  jumpreturn delay
         dm = 'nnn':  no broadband decoupling of heteronucleus during
                      acquisition; 
              'nny':  broadband heteronuclear decoupling during acquisition

      phase = 1,2:  hypercomplex experiment with F1 quadrature (complex F1-FT)
  */



#include <standard.h>


pulsesequence()
{
/* VARIABLE DECLARATION */
   double          ss,at,
                   pwx2lvl,
                   pwx2,
                   j,sw1,
                   bird,
		   total_time,
		   decup_time,
		   tau;
   int             iphase,
                   decup_ok;
    loadtable("11hmqc3rf");

/* LOAD VARIABLES */
   pwx2lvl = getval("pwx2lvl");
   pwx2 = getval("pwx2"); at=getval("at");
   sw1=getval("sw1");
   j = getval("j");
   ss = getval("ss");
   iphase = (int) (getval("phase") + 0.5);
   tau = getval("tau");


/* INITIALIZE VARIABLES */
   if (j > 0.0)
   {
      bird = 1.0 / (2.0*j);
   }
   else
   {
      bird = 0.0;
   }

/* Check for correct DM settings */
   if ((dm[A] == 'y') || (dm[B] == 'y') )
   {
      printf("DM must be set to either 'nny' or 'nnn'.\n");
      abort(1);
   }

/* Check for ROF1 minimum value */
   if ((rof1 < 9.9e-6) && (ix == 1))
      printf("Warning:  ROF1 is less than 10 us\n");

/* Check for correct decoupler duty cycle */
      decup_time = 2.0*pwx2;
   if (dm[C] == 'y')
      decup_time = decup_time + at;

   total_time =  d1 + at + 2.0*bird ;

   decup_ok = (((decup_time / total_time) < 0.2) && (dpwr2 < 50.0));
   if (!decup_ok)
   {
      printf("Decoupler duty cycle must be under 20%.\n");
      abort(1);
   }


/* DETERMINE STEADY-STATE MODE */
   if (ss < 0)
   {
      ss = (-1) * ss;
   }
   else
   {
      if ((ss > 0) && (ix == 1))
      {
	 ss = ss;
      }
      else
      {
	 ss = 0;
      }
   }
   initval(ss, ssval);
   initval(ss, ssctr);
   getelem(t3,ct,v3);
   if (iphase == 2)
      incr(v3);
   getelem(t7,ct,oph);

   if ((iphase == 1) || (iphase == 2))
   {
      initval(2.0*(double)((int)(d2*sw1+0.5)%2),v14);
      add(v3,v14,v3); add(oph,v14,oph);
   }

/* BEGIN ACTUAL PULSE SEQUENCE CODE */
   status(A);
      rlpower(pwx2lvl, DO2DEV);
      rlpower(tpwr,TODEV);
      hsdelay(d1);
   status(B);
      rcvroff();
      rgpulse(pw, t1, rof1, rof1);
      delay(tau - 2.0*rof1);
      rgpulse(pw-p1, t2, rof1, rof1);
      txphase(t4);
      dec2phase(v3);
      delay(bird - 0.5*pw - 0.5*pwx2 - 2.0*rof1);
         dec2rgpulse(pwx2, v3, rof1, 0.0);
         if (d2>0) delay(d2/2.0 - pw -tau -rof1  - 2.0*pwx2/3.1414);
         dec2phase(t6);
	 rgpulse(pw, t4, rof1, rof1);
         delay(2.0*tau -2.0*rof1);
	 rgpulse(pw, t5, rof1, rof1);
         if (d2>0) delay(d2/2.0 - pw -tau -2.0*rof1 - 2.0*pwx2/3.1414);
	 dec2rgpulse(pwx2, t6, rof1, 0.0);
      delay(rof2);
      rlpower(dpwr2, DO2DEV);
      delay(bird - 4.2e-6 - 0.5*pwx2 - rof2);
      rcvron();
   status(C);
}
/*
t1 = 0
t2 = 2
t3 = 0 2
t4 = 0
t5 = 2
t6 = 0 0 2 2
t7 = 0 2 2 0
*/

