/* noesyskbx -  2D cross relaxation experiment using sklbax "read" pulse
 
   (includes : 1. FAD-hypercomplex
               2. correction for d2 > 0          )

       mix = mixing time.
     phase = 1,2: gives HYPERCOMPLEX phase-sensitive experiment;
             3 :  gives TPPI phase-sensitive experiment;
    p1     = soft pulse width entered as delay in us - calibrate with
              the first increment
             (*** p1 calibrated using the first increment with this
              sequence is generally less than that calibrated using
              sklbax sequence ***)
    softpwr = power for soft pulse
    trim   = trim pulse follows the read pulse
  phaseinc = phase increment (in degree) to compensate phase difference
             between the soft pulse (p1) and the hard pulse (pw)


          -------    minimum 16 scans recommended    ------

                                Sklenar and Bax, JMR, 74,469(1987)

   contact : Krish Krishnamurthy, Florham Park, Sept. 1990               */


#include <standard.h>

pulsesequence()
{
double arraydim,
       ni,
       softpwr,
       trim,
       phaseinc,
       mix,
       phase;

int    iphase;

char   sspul[MAXSTR];

   stepsize(0.5,TODEV);
   phaseinc = getval("phaseinc");
   phaseinc = 2.0*phaseinc;
    if (phaseinc < 0.0) phaseinc=720+phaseinc;
   initval(phaseinc, v6);

   ni = getval("ni");
   arraydim = getval("arraydim");
   trim = getval("trim");
   softpwr = getval("softpwr");
    initval(softpwr, v10);
   tpwr = getval("tpwr");
    initval(tpwr, v11);
   mix = getval("mix");
   phase = getval("phase");
   getstr("sspul",sspul);
   iphase = (int)(phase + 0.5);

   loadtable("noesyskbx");

   getelem(t1,ct,v1);
   getelem(t6,ct,oph);
   
   if (iphase == 2)
      incr(v1);      /* hypercomplex phase shift */

   if (iphase == 3)                                  /*TPPI */
     initval((double)((int)((ix-1)/(arraydim/ni)+1.0e-6)), v8);
   else
   initval(2.0*(double)(((int)(d2*getval("sw1")+0.5)%2)),v8);
                                             /*FAD - Hypercomplex */
   add(v1,v8,v1);   /*TPPI or FAD-Hypercomplex */

   if ((iphase == 1) || (iphase == 2))
    add(oph,v8,oph);      /*FAD - Hypercomplex */

 /*ACTUAL PULSE SEQUENCE BEGINS */

   status(A);
   power(zero,DODEV);          /*reduced decoupler power */
   if (sspul[0] == 'y')
    {
     rgpulse(200*pw, zero, rof1,0.0);
     rgpulse(200*pw, one, 0.0,rof2);
    }                                             
     hsdelay(d1);

   status(B);
      rgpulse(pw, v1, rof1, 1.0e-6);
      if (d2 > 0)
       delay(d2 - rof1 - 1.0e-6 - (4*pw/3.1416));     /*corrected evolution */
      rgpulse(pw, t2, rof1, 1.0e-6);

   status(C);
      power(v10,TODEV);
      xmtrphase(v6);
      hsdelay(mix);

   status(D);                                     
      rgpulse(p1,t5,4.0e-6,0.0);
      power(v11,TODEV);
      xmtrphase(zero);
      rgpulse(pw,t3,rof1,0.0);
      rgpulse(trim,t4,2e-6,rof2);
}

/*   noesyskbx phase table:

     t1 = 0                                               - first pulse
     t2 = [0 2]4                                          - second pulse
     t3 = (0 2)4 (1 3)4                                   - hard pulse
     t4 = (1 1 3 3)2 (2 2 0 0)2                           - trim pulse
     t5 = [3 0 1 2]8                                      - soft pulse
     t6 = (0 2)2 (2 0)2 (1 3)2 (3 1)2                     - receiver */
