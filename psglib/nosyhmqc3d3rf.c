/* nosyhmqc3d3rf -  NOESY-HMQC 3D sequence 
                - written in hypercomplex phase sensitive mode only
                  uses second decoupler for X 
  Sequence:

  NOESY-HMQC:

  status :    A-|-B|----------------------C-----------------|-D
     1H  :    90-t1-90-mix-90-1/2J-        180        -1/2J- Acq (t3)
      X  :       BB                90-t2/2-   -t2/2-90       BB
 phtable :    t1    t3     t3      t2       t3      t3       t4 or t5

  Parameters:

        d2 = First evolution time
        d3 = Second evolution time
       mix = NOESY mixing time.
    pwx2lvl = power level for X pulses
    pwx2    = 90 degree X pulse
    jxh    = X-H coupling constant
    dpwr2   = power level for X decoupling 
    tpwr   = power level for H pulses
    pw     = 90 degree H pulse
     phase =   1,2: gives HYPERCOMPLEX (t1) acquisition;
     ni    = number of t1 increments
     phase2=   1,2: gives HYPERCOMPLEX (t2) acquisition;
     ni2   = number of t2 increments
    satmode  = flag to do xmtr presaturation at satfrq and satpwr
              during relaxation (satdly) and mixing (mix) periods
    satfrq  = saturation frequency for xmtr presaturation
    satpwr  = saturation power for xmtr presaturation
    satdly  = saturation period follows D1
    sspul = 'y': selects for trim(x)-trim(y) sequence at start of pulse sequence
    nosyflg= 'n': turns off NOESY part of the sequence 
    hmqcflg= 'n': turns off HMQC part of the sequence 

  Plane processing:

     array = phase2,phase

        np x ni  :  wft2d('ni',#i,1,0,0,0,0,0,0,0,0,0,-1,0,0,0,0,0)
        np x ni2 : wft2d('ni2',#i,1,0,0,0,0,0,0,0,0,0,0,0,-1,0,0,0)

     array = phase,phase2

        np xni   :  wft2d('ni',#i,1,0,0,0,0,0,0,0,0,0,0,0,-1,0,0,0)
        np x ni2 : wft2d('ni2',#i,1,0,0,0,0,0,0,0,0,0,-1,0,0,0,0,0)

  phase table:

        t1 = 0 2 0 2 1 3 1 3  - first proton pulse
        t2 = 0 0 2 2 1 1 3 3  - first carbon pulse
        t3 = 0 0 0 0 1 1 1 1  - all other pulses
        t4 = 0 2 2 0 1 3 3 1  - receiver
        t5 = 0 2 0 2 1 3 1 3  - receiver if hmqcflg='n'

              Feb 07 1991    (vvk)
    revised:  March 01, 1991 (vvk)*/


#include <standard.h>

pulsesequence()
{
   double          satpwr,
                   satdly,
                   satfrq,
                   mix,
                   jxh,
                   pwx2lvl,
                   pwx2,
                   tau;
   int             phase2,
                   phase;
   char            sspul[MAXSTR],
                   satmode[MAXSTR],
                   nosyflg[MAXSTR],
                   hmqcflg[MAXSTR];

/* LOAD VARIABLES */
   satdly = getval("satdly");
   satfrq = getval("satfrq");
   satpwr = getval("satpwr");
   mix = getval("mix");
   pwx2lvl = getval("pwx2lvl");
   pwx2 = getval("pwx2");
   jxh = getval("jxh");
   tau = 1/(2.0*jxh);
   phase = (int) (getval("phase") + 0.5);
   phase2 = (int) (getval("phase2") + 0.5);
   getstr("sspul", sspul);
   getstr("satmode",satmode);
   getstr("nosyflg",nosyflg);
   getstr("hmqcflg",hmqcflg);


   initval(pwx2lvl,v6);
   initval(dpwr2,v7);
   initval(tpwr,v8);
   initval(satpwr,v9); 
  

/* CHECK CONDITIONS */
  

/* DETERMINE STEADY-STATE MODE */


/* STEADY-STATE PHASECYCLING

/* PHASECYCLE */

   loadtable("nosyhmqc3d3rf");    /* t1 = 0 2 0 2 1 3 1 3
                                  t2 = 0 0 2 2 1 1 3 3
                                  t3 = 0 0 0 0 1 1 1 1
                                  t4 = 0 2 2 0 1 3 3 1
                                  t5 = 0 2 0 2 1 3 1 3     */
   
   getelem(t1,ct,v1);
   getelem(t2,ct,v2);

   if (hmqcflg[A] == 'n')
    getelem(t5,ct,oph);
   else
    getelem(t4,ct,oph);

   initval(2.0*(double)(((int)(d2*getval("sw1")+0.5)%2)),v13);
   initval(2.0*(double)(((int)(d3*getval("sw2")+0.5)%2)),v14);

     if (phase == 2)
        incr(v1);
     if (phase2 == 2)
        incr(v2);
     add(v1,v13,v1);
     add(v2,v14,v2);

   add(oph,v13,oph);
   add(oph,v14,oph);

    

/* BEGIN THE ACTUAL PULSE SEQUENCE */
   status(A);              /* relaxation status */
   power(v8,TODEV);
   power(v7,DO2DEV);
   dec2phase(zero);
   if (sspul[A] == 'y')
   {
      rgpulse(200*pw, zero, rof1,0.0);
      rgpulse(200*pw, one, 0.0,rof2);
   }
      hsdelay(d1);

   if (satmode[A] == 'y')
     {
        if (satfrq != tof)
         offset(satfrq,TODEV);
        power(v9,TODEV);
        rgpulse(satdly,zero,rof1,rof1);
        power(v8,TODEV);
        if (satfrq != tof)
         {
          offset(tof,TODEV);
          delay(40.0e-6);
         }
     }
   else
     delay(satdly);
   
     if (nosyflg[A] != 'n')
      {
       rgpulse(pw, v1, rof1, 1.0e-6);

   status(B);
   if (satmode[B] == 'y')
     {
      if (d2>0)
       {
        if (satfrq != tof)
         offset(satfrq,TODEV);
        power(v9,TODEV);
        rgpulse(d2 - 10.2e-6 -3*rof1 - (4*pw/3.1416),zero,rof1,rof1);
        power(v8,TODEV);
        if (satfrq != tof)
         {
          offset(tof,TODEV);
          delay(40.0e-6);
         }
       }
     }
   else
     {
         if (d2>0)
          delay(d2 - rof1 - 1.0e-6 - (4*pw/3.1416)); 
     }
    status(C);
       rgpulse(pw, t3, rof1, 1.0e-6);
       if (satmode[C] == 'y')
        {
          if (satfrq != tof)
           offset(satfrq,TODEV);
          power(v9,TODEV);
          rgpulse(mix,zero,2.0e-6,rof1);
          power(v8,TODEV);
          if (satfrq != tof)
           {
            offset(tof,TODEV);
            delay(40.0e-6);
           }
         }
        else
         hsdelay(mix);
         power(v6,DO2DEV);
        rgpulse(pw, t3, rof1, 0.0);
       }
      else
       rgpulse(pw, v1, rof1, 0.0);


      if (hmqcflg[A] != 'n')
       {
        delay(tau - 0.5*pw - 0.5*pwx2 - rof1);
        txphase(t3);
        dec2rgpulse(pwx2,v2,rof1,0.0);
         dec2phase(t3);
         delay(d3/2.0-pw-(2.0*pwx2)/3.1414);
        rgpulse(2.0*pw, t3, 0.0, 0.0);
         dec2phase(t3);
         delay(d3/2.0-pw-(2.0*pwx2)/3.1414);
        dec2rgpulse(pwx2,t3,0.0,rof1);
        delay(tau - 4.2e-6 - 0.5*pwx2 - rof1);
       }
      power(v7,DO2DEV);
      delay(rof2);  

    status(D);                                    /* acquisition status */
}
