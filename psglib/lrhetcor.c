#ifndef LINT
static char     SCCSid[] = "@(#)lrhetcor.c";
#endif

/* lrhetcor - long-range heteronuclear chemical shift correlation with 
              BIRD pulse during refocusing delay and a 1-step
              low-pass J filter

  include:  (1) 1-step low-pass J filter
            (2) coupled or decoupled spectrum (for coupled spectrum
                the refocusing delay is automatically eliminated)
            (3) BIRD pulse during refocusing (for decoupled spectrum)


 Parameters:

      pw  -  90 degree pulse on the observe nucleus
    tpwr  -  transmitter power level; only for systems with a linear
             amplifier on the transmitter channel
      pp  -  proton 90 degree pulse on the decoupler channel
   pplvl  -  decoupler power level; only for systems with a linear
             amplifier on the decoupler channel; otherwise decoupler
             is turned to to full-power for pulses on systems that
             have bilevel decoupling capability
     dhp  -  decoupler power level during acquisition
    dpwr  -  decoupler power level during acquisition for systems with
	     linear amplifiers
    j1xh  -  one-bond heternuclear coupling constant
    jnxh  -  multiple-bond heteronuclear coupling constant
      nt  -  multiple of 8 (minimum required)
             multiple of 32 (minimum recommended)


 References:

  Gary Martin, Magn. Reson. Chem. 26, 28 (1988) 
  V.V. Krishnamurthy, J. Magn. Reson. 80, 280 (1988) 

  revised Oct 1991  (vvk)                                  */



#include <standard.h>

pulsesequence()
{
   double	j1xh,
         	jnxh,
         	pp,
	 	pplvl,
         	dly3,
         	dly4,
         	dly5;

/* Get new variables from parameter table */
   pp = getval("pp");
   j1xh = getval("j1xh");
   jnxh = getval("jnxh");
   dly3 = getval("dly3");
   dly4 = getval("dly4");
   dly5 = getval("dly5");
   if (newdecamp)
   {
      pplvl = getval("pplvl");
      initval(pplvl, v10);
      initval(dpwr, v11);
   }

/* Calculate delays  -  allows direct entry of dly3 and dly4 delays - preferred
                        if set to zero calculated from jnxh
                        if jnxh is also zero optimized for 8 Hz    
                         
                        allows direct entry of dly5
                        if set to zero calculated from j1xh - preferred
                        if j1xh is also zero optimized for 150 Hz   */
   if (dly3 == 0)
      {
       if (jnxh == 0)
          dly3 = 0.0625;
       else
          dly3 = 1.0/(2.0*jnxh);
      }
   if (dly4 == 0)
      {
       if (jnxh == 0)
          dly4 = 0.031;
       else
          dly4 = 1.0/(4.0*jnxh);
      }
   if (dly5 == 0)
      {
       if (j1xh == 0)
          dly5 = 0.0033;
       else
          dly5 = 1.0/(2.0*j1xh);
      }


   loadtable("lrhetcor");
 
                                /*lrhetcor table
                                  t1 = 0 2
                                  t2 = [0 2]16
                                  t3 = [0 1 2 3]2
                                  t4 = [0 1]8
                                  t5 = [0 1 2 3]2 [2 3 0 1]2
                                  t6 = [0 1 2 3]2
                                               */

   getelem(t2,ct,v1);
   add(v1,one,v2);

   getelem(t4,ct,v3);
   add(v3,one,v4);


/* ACTUAL PULSE SEQUENCE BEGINS */
/* Relaxation delay */
   status(A);
      if (dm[0] == 'y')
      {
         fprintf(stdout, "decoupler must be set as dm=nny or nnn\n");
         abort(1);
      }
 
      if (newdecamp)
      {
         power(v10, DODEV);
      }
      else
      {
         declvlon();		/* sets maximum dhp value */
      }
      hsdelay(d1);


   status(B);
      rcvroff();
      delay(2.0e-5);		/* amp unblanking time */
      decpulse(pp, zero);
      delay(dly5 - rof1);
      rgpulse(pw,t1,rof1,rof2);
      delay(d2/2.0);
      rgpulse(pw, v1, rof1, 0.0);
      rgpulse(2*pw, v2, 0.0, 0.0);
      rgpulse(pw, v1, 0.0, rof1);
      delay(d2/2.0 + rof2);
      delay(dly3);
      simpulse(pw, pp, zero, t3, rof1, rof2);
      if (dm[2] == 'y')
      {
        delay(dly4/2.0 - rof2);
        decpulse(pp,zero);
        delay(dly5 - rof1 - pw - 1.5 * pp);
        rgpulse(pw,v3,rof1,0.0);
        simpulse(2.0*pw, 2.0*pp, v4, one, 0.0,0.0);
        rgpulse(pw,v3,0.0,rof1);
        delay(dly5 - rof1 - pw - 1.5 * pp);
        decpulse(pp,zero);
        rcvron();
        if (newdecamp)
        {
          power(v11, DODEV);
        }
        else
        {
          declvloff();
        }
        delay(dly4/2.0);
        setreceiver(t5);
      }
      else
      {
        rcvron();
        if (newdecamp)
        {
           power(v11, DODEV);
        }
        else
        {
           declvloff();
        }
        setreceiver(t6);
      }
 
/* Observe period */
   status(C);
}
