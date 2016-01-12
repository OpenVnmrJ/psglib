/* hmqctoxy3d3rf - HMQC-TOCSY 3D sequence
              - written in hypercomplex phase sensitive mode only

   sequence:

   HMQC-TOCSY:

   status   : A|--------------B-----------------|-----C-----|--D
     1H     :   90-1/2J-        180        -1/2J-t2-spinlock-Acq (t3)
      X     :           90-t1/2-   -t1/2-90     -----BB---------
  phtable   :   t1      t3      t6       t6           t2      t4 or t5    

   Parameters:

         d2 = First evolution time
         d3 = second evolution time
        mix = TOCSY mixing time
     pwx2lvl = power level for X pulses
     pwx2    = 90 degree X pulse
     jxh    = X-H coupling constant
     dpwr2   = power level for X decoupling
     tpwr   = power level for H pulses
     pw     = 90 degree H hard pulse
     p1lvl  = power level for spinlock
     p1     = 90 degree H pulse for mlev17
     trim   = trim pulse preceeding mlev17
     phase  = 1,2: gives HYPERCOMPLEX (t1) acquisition;
     ni     = number of t1 increments
     phase2 = 1,2: gives HYPERCOMPLEX (t2) acquisition;
     ni2    = number of t2 increments
     satmode = 'y':  presaturation during relaxation period (satdly) with xmtr
              'yy':  presaturation during both relaxation period (satdly)
                     and null period (null) with xmtr
     satfrq = presaturation frequency
     satdly = saturation time during the relaxation period
     satpwr = saturation power for all periods of presaturation with xmtr
         hs = 'yn': homospoil pulse (hst) during the d1 relaxation delay
              'yy': homospoil pulse during both d1 and the nulling period (null)
    wdwfctr = multiplication "window" factor of p1
       null = nulling time for protons not attached to observed heteronucleus
    hmqcflg = 'n': turns off HMQC part of the sequence
    toxyflg = 'n': turns off TOCSY part of the sequence

  Plane processing:

    array = phase2,phase

      np x ni   :  wft2d('ni',#i,1,0,0,0,0,0,0,0,0,0,-1,0,0,0,0,0)
      np x ni2  : wft2d('ni2',#i,1,0,0,0,0,0,0,0,0,0, 0,0,1,0,0,0)

     array = phase,phase2

      np x ni   :  wft2d('ni',#i,1,0,0,0,0,0,0,0,0,0,0,0,-1,0,0,0)
      np x ni2  : wft2d('ni2',#i,1,0,0,0,0,0,0,0,0,0,1,0, 0,0,0,0)

  Phase table:

    t1 = 0 0 1 1                      - first proton pulse
    t2 = 0 0 1 1 2 2 3 3              - spinlock 
    t3 = 0 2 1 3                      - first carbon pulse
    t6 = 0 0 1 1                      - proton 180 and second carbon pulse
    t4 = 0 2 1 3                      - receiver
    t5 = 0 0 1 1                      - receiver if hmqcflg='n'
    t7 = 1 1 2 2                      - composite X 180 during BIRD pulse
    t8 = 2 2 3 3                      - second H 90 during BIRD pulse

   march 25,1991  (vvk)   */

#include <standard.h>


mleva()
{
   double wdwfctr,window;
   wdwfctr=getval("wdwfctr");
   window = (wdwfctr*p1);
   rgpulse(p1,v2,0.0,window);
   rgpulse(2*p1,v4,0.0,window);
   rgpulse(p1,v2,0.0,0.0);
}

mlevb()
{
   double wdwfctr,window;
   wdwfctr=getval("wdwfctr");
   window = (wdwfctr*p1);
   rgpulse(p1,v5,0.0,window);
   rgpulse(2*p1,v12,0.0,window);
   rgpulse(p1,v5,0.0,0.0);
}


pulsesequence()
{

/* VARIABLE DECLARATION */
   double  pwx2lvl,
           phase,
           phase2,
           pwx2,
           jxh,
           bird,
           null,
           satfrq,
           satdly,
           satpwr,
           cycles,
           trim,
           wdwfctr,
           window,
           p1lvl,
           at,
           mix;
    int iphase,
        iphase2;
    char satmode[MAXSTR],
         hmqcflg[MAXSTR],
         toxyflg[MAXSTR];


/* LOAD VARIABLES */

    loadtable("hmqctoxy3d");            
                           /* hmqctoxy3d phase table
                              t1 = 0 0 1 1           - H1 90 
                              t6 = 0 0 1 1           - H1 180 and 2nd X 90
                              t7 = 1 1 2 2 
                              t8 = 2 2 3 3          
                              t2 = 0 0 1 1 2 2 3 3   - spin lock
                              t3 = 0 2 1 3           - first X 90
                              t4 = 0 2 1 3           - receiver
                              t5 = 0 0 1 1           - receiver if hmqcflg=n
                           */
    satfrq = getval("satfrq");
    satdly = getval("satdly");
    satpwr = getval("satpwr");
    pwx2lvl = getval("pwx2lvl");
    mix = getval("mix");
    p1lvl = getval("p1lvl");
    trim = getval("trim");
    wdwfctr = getval("wdwfctr");
    window = (wdwfctr*p1);
    phase = getval("phase");
    phase2 = getval("phase2");
    pwx2 = getval("pwx2");
    null = getval("null");
    jxh = getval("jxh");
    at = getval("at");
    getstr("satmode",satmode);
    getstr("hmqcflg",hmqcflg);
    getstr("toxyflg",toxyflg);


/* INITIALIZE VARIABLES */
    if (jxh > 0.0)
       bird = 1.0/(2.0*jxh);
    else
       bird = 0.0;
 
    iphase = (int) (phase + 0.5);
    iphase2 = (int)(phase2 + 0.5);

       initval(2.0*(double)(((int)(d2*getval("sw1")+0.5)%2)),v14);
       initval(2.0*(double)(((int)(d3*getval("sw2")+0.5)%2)),v13);


    initval(pwx2lvl,v10);
    initval(dpwr2,v8);
    initval(p1lvl,v9);
    initval(satpwr,v7);
    initval(tpwr,v6);

/* SAFETY CHECKS    */

   if (satmode[A] == 'y')
    {
     if (((at+mix+trim)/(at+mix+trim+d1+satdly+null)) >= 0.25)
     {
      printf("Duty cycle too high.   Please check. \n");
      abort(1);
     }
    }
   else
    {
     if (((at+mix+trim)/(at+mix+trim+d1+null)) >= 0.25)
     {
      printf("Duty cycle too high.   Please check. \n"); 
      abort(1); 
     }
    }

   if ((at+mix+trim) >= 0.250)
     {
      printf("Decoupler is turned on for more than 250 ms. Please check.\n");
      abort(1);  
     }

   if ((dm2[A] == 'y') || (dm2[B] == 'y'))
     {
      printf("dm2 should be 'nnn' or 'nny'. \ n");
      abort(1);   
     } 

/* DETERMINE STEADY-STATE MODE */


/* STEADY-STATE PHASECYCLE */

    getelem(t1,ct,v1);
    getelem(t3,ct,v3);

   if (hmqcflg[0] == 'y')
    getelem(t4,ct,oph);
   else
    getelem(t5,ct,oph);


    if (iphase == 2)
       incr(v3);
    if (iphase2 == 2)
       incr(v1);

    add(v1,v13,v1);
    add(v3,v14,v3);

    add(oph,v14,oph);
    add(oph,v13,oph);
   

/* BEGIN ACTUAL PULSE SEQUENCE CODE */
    status(A);
       delay(1.0e-5); power(v10,DO2DEV);
       power(v6,TODEV);
       hsdelay(d1);

/* selective saturation period */
       if (satmode[A] == 'y')
         {
          offset(satfrq,TODEV);
          power(v7,TODEV);
          rgpulse(satdly,zero,4.0e-5,1.0e-5);
          offset(tof,TODEV);
          power(v6,TODEV); 
          delay(5.0e-5);
         }

/* if null is 0 eliminate bird inversion pulse */
/* bird pulse - 180 for protons not coupled to 13C */

    status(B);
    if (hmqcflg[0] == 'y')
     {
       if (null != 0.0)
         {
          rcvroff();
          rgpulse(pw,t6,1.0e-5,0.0);
          if (pw > pwx2)
            delay(bird - rof1 - 1.0e-6 - 1.5*pw - pwx2);
          else
            delay(bird - rof1 - 1.0e-6 - 2*pwx2 - 0.5*pw);
          dec2rgpulse(pwx2,t7,rof1,0.0);
          sim3pulse(2.0*pw,0.0,2.0*pwx2,t6,zero,t6,1.0e-6,0.0);
          dec2phase(t7);
          dec2rgpulse(pwx2,t7,1.0e-6,0.0);
          if (pw > pwx2)
            delay(bird - rof1 - 1.0e-6 - 1.5*pwx2 - pwx2);
          else
            delay(bird - rof1 - 1.0e-6 - 2*pwx2 - 0.5*pw); 
          rgpulse(pw,t8,rof1,rof2);
          rcvron();

/* nulling time for protons not coupled to 13C */

          if (satmode[1] == 'y')
            {
             offset(satfrq,TODEV);
             power(v7,TODEV);
             rgpulse(null,zero,4.0e-5,1.0e-5);
             offset(tof,TODEV);
             power(v6,TODEV);
             delay(4.0e-5);
            }
          else
             hsdelay(null);
         }

       rcvroff();
       rgpulse(pw,v1,rof1,0.0);
       delay(bird - 0.5*pw - 0.5*pwx2 - rof1);
       txphase(t6);
       if (d2 <= 2*pw)
       {
          dec2phase(v3);
          delay(rof1);
          if (pwx2 >= (pw - 0.5e-6))
          {
            dec2on();
            delay(pwx2 - (2*pw - d2)/2);
            xmtron();
            delay(pw - d2/2);
            dec2off();
            dec2phase(t6);
            delay(d2);
            dec2on();
            delay(pw - d2/2);
            xmtroff();
            delay(pwx2 - (2*pw - d2)/2);
            dec2off();
          }
          else
          {
            xmtron();
            delay((2*pw - d2)/2 - pwx2);
            dec2on();
            delay(pwx2);
            dec2off();
            dec2phase(t6);
            delay(d2);
            dec2on();
            delay(pwx2);
            dec2off();
            delay((2*pw - d2)/2 - pwx2);
            xmtroff();
          }
        }
        else
        {
          dec2rgpulse(pwx2,v3,rof1,0.0);
          delay(d2/2.0 - pw + 0.5e-6 -(2.0*pwx2/3.1414));
          dec2phase(t6);
          rgpulse(2.0*pw,t6,0.0,0.0);
          delay(d2/2.0 - pw + 0.5e-6 -(2.0*pwx2/3.1414));
          dec2rgpulse(pwx2,t6,0.0,0.0);
        }
  
       delay(bird);
      }
     else
      {
       rcvroff();
       rgpulse(pw,v1,rof1,0.0);
      }
       power(v8,DO2DEV);
       power(v9,TODEV);
       dec2phase(zero);          /* to allow x decoupling */

   status(C);
     if (toxyflg[0] == 'y')
      {
       getelem(t2,ct,v2);
       add(one,v2,v4);
       add(one,v4,v5);
       add(one,v5,v12);
       txphase(v4);
        if (d3>0.0)
         {
          if (hmqcflg[0]=='y')
          delay(d3 - 9.4e-6);
          else
          delay(d3 - 9.4e-6 - (2*pw/3.1416));
         }
        else
          delay(d3);
       

/* calculate and initialize loop counter */
       cycles = (mix-trim)/(64.66*p1+32*window);
       cycles = 2.0*(double)(int)(cycles/2.0);
       initval(cycles,v11);
       if (cycles > 1.0)
        {
         rgpulse(trim,v4,5.0e-6,0.0);
          starthardloop(v11);
             mleva(); mlevb(); mlevb(); mleva();
             mlevb(); mlevb(); mleva(); mleva();
             mlevb(); mleva(); mleva(); mlevb();
             mleva(); mleva(); mlevb(); mlevb();
             rgpulse(0.66*p1,v4,0.0,0.0);
          endhardloop();
        }
      }
       delay(rof2);
       rcvron();
       txphase(zero);
    status(D);
}
