/* nosytoxy3d -  NOESY-TOCSY 3D sequence 
                 - written in hypercomplex phase sensitive mode only
  
  Sequence:

  NOESY-TOCSY:

  status : A-|------B--------|----C------|--D
     1H  :    90-t1-90-mix-90-t2-spinlock-Acq (t3)
                                 (MLEV17)
  phtable:    t1    t2     t3       t4     t5 or t3

  Parameters:

        d2 = First evolution time
        d3 = Second evolution time
       mix = NOESY mixing time.
     slmix = spin-lock mixing time
      trim = trim pulse before mlev17
     trim2 = trim pulse after mlev17 (optional)
     p1lvl = power level for mlev17
     p1    = 90 degree pulse for mlev17
      tpwr = power level for hard pulses
        pw = 90 degree hard pulse
     phase =   1,2: gives HYPERCOMPLEX (t1) acquisition;
     ni    = number of t1 increments
     phase2=   1,2: gives HYPERCOMPLEX (t2) acquisition;
     ni2   = number of t2 increments
   satflg  = 'yy' uses obs xmtr for presat at satfrq and satpwr
              during satdly and mix
   satfrq  = saturation frequency for xmtr presaturation
   satpwr  = saturation power for xmtr presaturation
   satdly  = saturation period follows D1
     sspul = 'y': selects for HS-90-HS sequence at start of pulse sequence
   nosyflg = 'n' - turns off NOESY part of the sequence
   toxyflg = 'n' - turns off TOCSY part of the sequence

  Plane processing:

     array = phase2,phase
          np x ni   :  wft2d('ni', #i,1,0,0,0,0,0,0,0,0,0,-1,0,0,0,0,0)
          np x ni2  :  wft2d('ni2',#i,1,0,0,0,0,0,0,0,0,0,0,0, 1,0,0,0)

     array = phase,phase2
          np x ni   :  wft2d('ni', #i,1,0,0,0,0,0,0,0,0,0,0,0,-1,0,0,0)
          np x ni2  :  wft2d('ni2',#i,1,0,0,0,0,0,0,0,0,0,1,0, 0,0,0,0)

  Phase table:

       t1 = (0 2)4                - first 90
       t2 = 0                     - second 90
       t3 = 0 2 1 3 2 0 3 1       - third 90 (and receiver if toxyflg='n')
       t4 = 0 2 1 3 2 0 3 1       - MLEV-17
       t5 = 0 0 1 1 2 2 3 3       - receiver


   3D processing:

     array = phase2, phase
   
       1  0  0  0  0  0  0  0
       0  0 -1  0  0  0  0  0
       0  0  0  0 -1  0  0  0
       0  0  0  0  0  0  1  0
       1  0  0  0  0  0 -1  0

     array = phase, phase2

       1  0  0  0  0  0  0  0 
       0  0  0  0 -1  0  0  0
       0  0 -1  0  0  0  0  0
       0  0  0  0  0  0  1  0
       1  0  0  0  0  0 -1  0
           March 15, 1991 (vvk)
  revised: March 18, 1991 (vvk)  */





#include <standard.h>

mleva()
{
  double wdwfctr,window;
  wdwfctr=getval("wdwfctr");
  window=(wdwfctr*p1);
  txphase(v4); delay(p1);
  xmtroff(); delay(window); xmtron();
  txphase(v5); delay(2*p1);
  xmtroff(); delay(window); xmtron();
  txphase(v4); delay(p1);
}

mlevb()
{
  double wdwfctr,window; 
  wdwfctr=getval("wdwfctr"); 
  window=(wdwfctr*p1); 
  txphase(v6); delay(p1);
  xmtroff(); delay(window); xmtron();
  txphase(v7); delay(2*p1);
  xmtroff(); delay(window); xmtron();
  txphase(v6); delay(p1);
} 
 
pulsesequence()
{
   double          p1lvl,
                   satpwr,
                   satdly,
                   satfrq,
                   mix,
                   slmix,
                   trim,
                   trim2,
                   window,
                   cycles;
   int             iphase2,
                   iphase;
   char            sspul[MAXSTR],
                   satflg[MAXSTR],
                   nosyflg[MAXSTR],
                   toxyflg[MAXSTR];


/* LOAD VARIABLES */
   satdly = getval("satdly");
   satfrq = getval("satfrq");
   satpwr = getval("satpwr");
   mix = getval("mix");
   p1lvl=getval("p1lvl");
   slmix = getval("slmix");
   trim = getval("trim");
   trim2=getval("trim2");
   window=(getval("wdwfctr")*p1);
   iphase = (int) (getval("phase") + 0.5);
   iphase2 = (int) (getval("phase2") + 0.5);
   getstr("sspul", sspul);
   getstr("satflg",satflg);
   getstr("nosyflg",nosyflg);
   getstr("toxyflg",toxyflg);


   initval(tpwr,v8);
   initval(satpwr,v9); 
   initval(p1lvl,v10);


/* CHECK CONDITIONS */
  

/* DETERMINE STEADY-STATE MODE */


/* STEADY-STATE PHASECYCLING */

/* PHASECYCLE */

   loadtable("nosytoxy3d");         /*nosytoxy3d phase table
                                      t1 = (0 2)4
                                      t2 = 0
                                      t3 = 0 2 1 3 2 0 3 1
                                      t4 = 0 2 1 3 2 0 3 1
                                      t5 = 0 0 1 1 2 2 3 3
                                     */
   
   getelem(t1,ct,v1);
   getelem(t3,ct,v3);

   if (nosyflg[A] == 'n')
    getelem(t3,ct,oph);
   else
    getelem(t5,ct,oph);

   initval(2.0*(double)(((int)(d2*getval("sw1")+0.5)%2)),v13);
   initval(2.0*(double)(((int)(d3*getval("sw2")+0.5)%2)),v14);

     if (iphase == 2)
        incr(v1);
     if (iphase2 == 2)
        decr(v3);
     add(v1,v13,v1);
     add(v3,v14,v3);

   add(oph,v13,oph);
   add(oph,v14,oph);

    

/* BEGIN THE ACTUAL PULSE SEQUENCE */
   status(A);              /* relaxation status */
   power(v8,TODEV);
   if (sspul[A] == 'y')
   {
      hsdelay(hst + 0.001);
      rgpulse(pw, zero, rof1,rof2);
      hsdelay(hst+d1);
   }
   else
      hsdelay(d1);

   if (satflg[A] == 'y')
     {
        if (satfrq != tof)
         offset(satfrq,TODEV);
        power(v9,TODEV);
        rgpulse(satdly,zero,rof1,rof1);
        power(v8,TODEV);
        if (satfrq != tof)
          offset(tof,TODEV);
        delay(40.0e-6);
     }
   else
     delay(satdly);
   
   status(B);
     if (nosyflg[A] != 'n')
      {
       rgpulse(pw, v1, rof1, 1.0e-6);

         if (d2>0)
          delay(d2 - rof1 - 1.0e-6 - (4*pw/3.1416)); 

       rgpulse(pw, t2, rof1, 1.0e-6);
       if (satflg[B] == 'y')
        {
          if (satfrq != tof)
           offset(satfrq,TODEV);
          power(v9,TODEV);
          rgpulse(mix,zero,2.0e-6,rof1);
          power(v8,TODEV);
          if (satfrq != tof)
           offset(tof,TODEV);
          delay(40.0e-6);
         }
        else
         hsdelay(mix);
       }

       rgpulse(pw, v3, rof1, 1.0e-6);

    status(C);
      if (toxyflg[A] != 'n')
       {
        power(v10,TODEV);
        getelem(t4,ct,v4);
        add(one,v4,v5);
        add(one,v5,v6);
        add(one,v6,v7);
        rcvroff();
        txphase(v5);
        if (d3 > 0.0)
         delay(d3 - 3.2e-6 - (2*pw/3.14159));
        else
         delay(d3);

        cycles = (slmix-trim-trim2)/(64.66*p1+32*window);
        cycles = 2.0*(double)(int)(cycles/2.0);
        initval(cycles,v2);

        xmtron();
        delay(trim);
        if (cycles > 1.0)
         {
          starthardloop(v2);
            mleva(); mleva(); mlevb(); mlevb();
            mlevb(); mleva(); mleva(); mlevb();
            mlevb(); mlevb(); mleva(); mleva();
            mleva(); mlevb(); mlevb(); mleva();
            txphase(v5);
            delay(0.66*p1);
          endhardloop();
         }
         if (trim2 != 0.0)
          {
           txphase(v5);
           delay(trim2);
         }
         xmtroff();

         delay(rof2);
         rcvron();
        }
        else
         delay(rof2-1.0e-6);

    status(D);                                    /* acquisition status */
}
