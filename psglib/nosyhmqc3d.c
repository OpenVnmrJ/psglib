/* nosyhmqc3d -  NOESY-HMQC 3D sequence 
                  with slp presaturation option
                - written in hypercomplex phase sensitive mode only

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
    pwxlvl = power level for X pulses
    pwx    = 90 degree X pulse
    jxh    = X-H coupling constant
    dpwr   = power level for X decoupling 
    tpwr   = power level for H pulses
    pw     = 90 degree H pulse
     phase =   1,2: gives HYPERCOMPLEX (t1) acquisition;
     ni    = number of t1 increments
     phase2=   1,2: gives HYPERCOMPLEX (t2) acquisition;
     ni2   = number of t2 increments
    satflg  = flag to do xmtr presaturation at satfrq and satpwr
              during relaxation (satdly) and mixing (mix) periods
    satfrq  = saturation frequency for xmtr presaturation
    satpwr  = saturation power for xmtr presaturation
    satdly  = saturation period follows D1
    slpflg = 'y' : slp presaturation
    h2off = water off-resonance frequency (tof - satfrq)
    satshape = slpsatd (tof is downfield to water)
               slpsatu (tof is upfield to water)
    sspul = 'y': selects for HS-90-HS sequence at start of pulse sequence
    nosyflg= 'n': turns off NOESY part of the sequence 
    hmqcflg= 'n': turns off HMQC part of the sequence 

  Plane processing:

     array = phase2,phase

        np x ni  :  wft2d('ni',#i,1,0,0,0,0,0,0,0,0,0,-1,0,0,0,0,0)
        np x ni2 : wft2d('ni2',#i,1,0,0,0,0,0,0,0,0,0,0,0,-1,0,0,0)

     array = phase,phase2

        np xni   :  wft2d('ni',#i,1,0,0,0,0,0,0,0,0,0,0,0,-1,0,0,0)
        np x ni2 : wft2d('ni2',#i,1,0,0,0,0,0,0,0,0,0,-1,0,0,0,0,0)

  3D processing:

     array = phase2,phase

         1  0  0  0  0  0  0  0
         0  0 -1  0  0  0  0  0
         0  0  0  0  1  0  0  0
         0  0  0  0  0  0 -1  0
         1  0  0  0  0  0 -1  0

     array = phase, phase2

         1  0  0  0  0  0  0  0 
         0  0  0  0 -1  0  0  0 
         0  0  1  0  0  0  0  0 
         0  0  0  0  0  0 -1  0
         1  0  0  0  0  0 -1  0

  phase table:

        t1 = 0 2 0 2 1 3 1 3  - first proton pulse
        t2 = 0 0 2 2 1 1 3 3  - first carbon pulse
        t3 = 0 0 0 0 1 1 1 1  - all other pulses
        t4 = 0 2 2 0 1 3 3 1  - receiver
        t5 = 0 2 0 2 1 3 1 3  - receiver if hmqcflg='n'

              Feb 07 1991    (vvk)
    revised:  March 01, 1991 (vvk)
              Feb 10, 1992*/


#include <standard.h>

pulsesequence()
{
   double          satpwr,
                   cycle1,
                   cycle2,
                   shape_pw,
                   h2off,
                   satdly,
                   satfrq,
                   mix,
                   jxh,
                   pwxlvl,
                   pwx,
                   tau;
   int             phase2,
                   phase;
   char            sspul[MAXSTR],
                   satflg[MAXSTR],
                   satshape[MAXSTR],
                   slpflg[MAXSTR],
                   nosyflg[MAXSTR],
                   hmqcflg[MAXSTR];

/* LOAD VARIABLES */
   satdly = getval("satdly");
   satfrq = getval("satfrq");
   satpwr = getval("satpwr");
   h2off = getval("h2off");
   mix = getval("mix");
   pwxlvl = getval("pwxlvl");
   pwx = getval("pwx");
   jxh = getval("jxh");
   tau = 1/(2.0*jxh);
   phase = (int) (getval("phase") + 0.5);
   phase2 = (int) (getval("phase2") + 0.5);
   getstr("sspul", sspul);
   getstr("satflg",satflg);
   getstr("satshape",satshape);
   getstr("slpflg",slpflg);
   getstr("nosyflg",nosyflg);
   getstr("hmqcflg",hmqcflg);


   if ((slpflg[0] == 'y') && (h2off != 0.0))
   {
   shape_pw = 10/h2off;
   if (shape_pw < 0.0)
    shape_pw = -shape_pw;

   cycle1 = (satdly/(shape_pw + 15.4e-6));
   cycle1 = 2.0*(double)(int)(cycle1/2.0);
   initval(cycle1,v8);

   cycle2 = (mix/(shape_pw + 15.4e-6));
   cycle2 = 2.0*(double)(int)(cycle2/2.0);
   initval(cycle2,v9);
   }

/* CHECK HARDWARE CONFIGURATION */

   if (rfwg[0] != 'y')
   {
    if ((slpflg[0] == 'y') || (slpflg[2] == 'y'))
    {
     printf("NOTE: slpflg=y requires PPM in Channel 1. psg aborted.\n");
     abort(1);
    }
   }


/* PHASECYCLE */

   loadtable("nosyhmqc3d");    /* t1 = 0 2 0 2 1 3 1 3
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
   rlpower(tpwr,TODEV);
   rlpower(dpwr,DODEV);
   decphase(zero);
   if (sspul[A] == 'y')
   {
      rgpulse(200*pw, zero, rof1,0.0);
      rgpulse(200*pw, one, 0.0, rof1);
   }

   hsdelay(d1);

   if (satflg[0] == 'y')
     {
        rlpower(satpwr,TODEV);
        if (slpflg[0] == 'n')
        {
        if (satfrq != tof)
         offset(satfrq,TODEV);
        rgpulse(satdly,zero,rof1,rof1);
        if (satfrq != tof)
          offset(tof,TODEV);
        }
        if (slpflg[0] == 'y')
        {
         if (h2off != 0.0)
         {
         rcvroff();
         delay(rof1);
         starthardloop(v8);
          shaped_pulse(satshape,shape_pw,zero,0.0,1e-6);
         endhardloop();
         delay(rof2);
         rcvron();
         }
         else
         {
         if (satfrq != tof)
          offset(satfrq,TODEV);
         shaped_pulse(satshape,satdly,zero,rof1,rof1);
         if (satfrq != tof)
          offset(tof,TODEV);
         }
        }
        rlpower(tpwr,TODEV);
        delay(40.0e-6);
     }
   
  status(B);
     if (nosyflg[0] != 'n')
      {
       rgpulse(pw, v1, rof1, 1.0e-6);
       if (d2>0)
        delay(d2 - rof1 - 1.0e-6 - (4*pw/3.1416)); 
       else
        delay(d2);
       rgpulse(pw, t3, rof1, 1.0e-6);
  
  status(C);
       if (satflg[2] == 'y')
        {
          rlpower(satpwr,TODEV);
          if (slpflg[2] == 'n')
          {
          if (satfrq != tof)
           offset(satfrq,TODEV);
          rgpulse(mix,zero,2.0e-6,rof1);
          if (satfrq != tof)
           offset(tof,TODEV);
          }
          if (slpflg[2] == 'y')
          {
           if (h2off != 0.0)
           {
           rcvroff();
           delay(rof1);
           starthardloop(v9);
            shaped_pulse(satshape,shape_pw,zero,0.0,1e-6);
           endhardloop();
           delay(rof2);
           rcvron();
           }
           else
           {
            if (satfrq != tof)
             offset(satfrq,TODEV);
            shaped_pulse(satshape,mix,zero,2.0e-6,rof1);
            if (satfrq != tof)
             offset(tof,TODEV);
           }
          }
          rlpower(tpwr,TODEV);
          delay(40.0e-6);
         }
        else
         hsdelay(mix);
         rlpower(pwxlvl,DODEV);
        rgpulse(pw, t3, rof1, 0.0);
       }
      else
       rgpulse(pw, v1, rof1, 0.0);


      if (hmqcflg[A] != 'n')
       {
        rcvroff();
        delay(tau - 0.5*pw - 0.5*pwx - rof1);
        txphase(t3);
        decrgpulse(pwx,v2,rof1,0.0);
         decphase(t3);
         delay(d3/2.0);
        rgpulse(2.0*pw, t3, 0.0, 0.0);
         decphase(t3);
         delay(d3/2.0);
        decrgpulse(pwx,t3,0.0,rof1);
        rcvron();
        delay(tau - 4.2e-6 - 0.5*pwx - rof1);
       }
      rlpower(dpwr,DODEV);
      delay(rof2);  

    status(D);                                    /* acquisition status */
}
