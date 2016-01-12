/* shphmqc- heteronuclear multiple-quantum coherence with hypercomplex
         implementation, steady-state transients for either all or
         just the first t1 increment, pulse+receiver phasecycling during
         steady-state transients, proton pulses are shaped pulses
          REQUIRES WAVEFORM GENERATOR
          USES THE SECOND CHANNEL FOR X PULSES AND DECOUPLING

   ref: Summers et al., j. amer. chem. soc. 108:4285-4294 (1986)


   Parameters:

        dof = frequency for decoupling during acquisition
     pwxlvl = power level for decoupler pulses
        pwx = 90 degree decoupler pulse length for the heteronucleus
          j = one-bond heteronuclear coupling constant (in Hz)
         dm = 'nnn':  no broadband decoupling of heteronucleus during
                      acquisition; 
              'nny':  broadband heteronuclear decoupling during acquisition
      shape= shapelib file for shaped pulse
         tpwr= power level for shaped pulses
           pw= pulse width for both shaped pulses
      phase = 1,2:  hypercomplex experiment with F1 quadrature (complex F1-FT)
  */



#include <standard.h>


pulsesequence()
{
/* VARIABLE DECLARATION */
   double          ss,at,
                   pwxlvl,
                   pwx,
                   j,
                   bird,
		   total_time,
		   decup_time;
   int             iphase,
                   decup_ok;
   char            
                  shape[MAXSTR];
    loadtable("shphmqc");

/* LOAD VARIABLES */
   pwxlvl = getval("pwxlvl");
   pwx = getval("pwx"); at=getval("at");
   j = getval("j");
   ss = getval("ss"); getstr("shape",shape);
   iphase = (int) (getval("phase") + 0.5);


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
      decup_time = 2.0*pwx;
   if (dm[C] == 'y')
      decup_time = decup_time + at;

   total_time =  d1 + at + 2.0*bird ;

   decup_ok = (((decup_time / total_time) < 0.2) && (dpwr < 50.0));
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
      initval(2.0*(double)((int)(d2*getval("sw1")+0.5)%2),v14);
      add(v3,v14,v3); add(oph,v14,oph);
   }

/* BEGIN ACTUAL PULSE SEQUENCE CODE */
   status(A);
      rlpower(pwxlvl, DODEV);
   if (tpwr>50)
   {
      printf("TPWR should be lower for shaped pulse.\n");
      abort(1);
   }
      rlpower(tpwr,TODEV);
      hsdelay(d1);
   status(B);
      shaped_pulse(shape,pw,t1,rof1,rof1);
      rlpower(tpwr+6,TODEV);
      delay(bird -4.2e-6 - 0.5*pwx - rof1);
         decpulse(pwx, v3 );
         if (d2>0) delay(d2/2.0  -rof1 - pw/2  - 2.0*pwx/3.1414);
         decphase(v5);
      shaped_pulse(shape,pw,t4,rof1,rof1);
         if (d2>0) delay(d2/2.0 - pw/2 -rof1 - 2.0*pwx/3.1414);
	 decpulse(pwx, t6);
      delay(rof2);
      rlpower(dpwr, DODEV);
      delay(bird - 4.2e-6 - 0.5*pwx - rof1);
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

