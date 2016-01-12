
/* noesyhmqc3rf- noesy followed by hmqc using 3rd channel for X 
                 F1=1H, F2=1H
   Parameters:

     satmode = 'ynnnn': presaturation during relaxation period (satdly)
              'ynyyn': presaturation during both relaxation period (satdly),
                       d3,and mix periods 
     satfrq = presaturation frequency
     satdly = saturation time during the relaxation period
     satpwr = saturation power for all periods of presaturation with xmtr
         hs = 'yn':  homospoil pulse (hst) during the  relaxation delay
         pw = 90 degree xmtr pulse length for protons (the observed nucleus)
     pwx2lvl = power level for decoupler pulses
        pwx2 = 90 degree decoupler pulse length for the heteronucleus
          j = one-bond heteronuclear coupling constant (in Hz)
    contact - G.Gray (palo alto)
   revision - modified from hmqc.c
     */

#include <standard.h>

pulsesequence()
{
/* VARIABLE DECLARATION */
   double          
                   phase,
                   pwx2lvl,
                   pwx2,
                   j,
                   bird,
                   satfrq,
                   satdly,
                   satpwr,
                   mix;
   char            satmode[MAXSTR];


/* LOAD VARIABLES */
   satfrq = getval("satfrq");
   satdly = getval("satdly");
   satpwr = getval("satpwr");
   phase  = getval("phase");
   pwx2lvl = getval("pwx2lvl");
      mix = getval("mix");
   pwx2 = getval("pwx2");
   j = getval("j");
   getstr("satmode", satmode);
   loadtable("noesyhmqc");     /* read phase table */

/* INITIALIZE VARIABLES */
   if (j > 0.0)
   {
      bird = 1.0 / (2.0*j);
   }
   else
   {
      bird = 0.0;
   }

   initval(pwx2lvl, v10);
   initval(satpwr, v7);
   initval(tpwr, v13);
   initval(dpwr2, v6);


/* CHECK CONDITIONS */


/* Check for ROF1 minimum value */
   if ((rof1 < 9.9e-6) && (ix == 1))
      printf("Warning:  ROF1 is less than 10 us\n");

/* BEGIN ACTUAL PULSE SEQUENCE CODE */
   status(A);
     power(v10, DO2DEV);
     hsdelay(d1);

/* selective saturation period */
     if (satmode[A] == 'y')
      {
         offset(satfrq, TODEV);
         power(v7, TODEV);
         rgpulse(satdly, zero, 4.0e-5, 1.0e-5);
         offset(tof, TODEV);
         power(v13, TODEV);
         delay(4.0e-5);
      }
  status(B);
    power(v6,DO2DEV);
    getelem(t1,ct,v9);
    if (phase > 1.5) incr(v9);
    rgpulse(pw, v9, rof1, 0.0);
    delay(d2 - rof1 -4.2e-6 -4.0*pw/3.1414);
    power(v10,DO2DEV);
  status(C);
    rgpulse(pw,t5,rof1,rof2);
     if (satmode[C] == 'y')
      {
         offset(satfrq, TODEV);
         power(v7, TODEV);
         rgpulse(mix, zero, 4.0e-5, 1.0e-5);
         offset(tof, TODEV);
         power(v13, TODEV);
         delay(4.0e-5);
      }
    else delay(mix);
    rgpulse(pw,t6,rof1,rof2);
  status(D);
      delay(bird);
         dec2rgpulse(pwx2, t2, rof1, 0.0);
	 rgpulse(2.0*pw, t3, rof1, 0.0);
	 dec2rgpulse(pwx2, t4, 2.0e-6, rof2);
      power(v6, DO2DEV);
      delay(bird);
   status(E);
    setreceiver(t7);
}

/* hmqcnoesy - phase table for isotope-filtered NOESY experiment

t1 = [0 2]2
t2 = 0
t3 = [0 2]4
t4 = 0 2
t5 = 0
t6 = 0
t7 = (0 2 2 0)2
t8 = 0               */
