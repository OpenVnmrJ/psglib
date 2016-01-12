/* nosyxfil -  NOESY-X-filtered 2D sequence 
                X-filter parallel to F1
               SLP PRESATURATION OPTION
                - written in hypercomplex phase sensitive mode only

  Sequence:

  NOESY-X-filter:

  status :    A-|-B|------------C- --------------|-D
     1H  :    90-t1-90-mix-90-1/2J-   180   -1/2J- Acq (t2)
      X  :       BB                90-   -90       BB
 phtable :    t1    t5     t6      t2 t6  t3       t4 or t7

  Parameters:

       mix = NOESY mixing time.
    pwxlvl = power level for X pulses
    pwx    = 90 degree X pulse
    jxh    = X-H coupling constant
    dpwr   = power level for X decoupling 
    tpwr   = power level for H pulses
    pw     = 90 degree H pulse
    satflg  = flag to do xmtr presaturation at satfrq and satpwr
              during relaxation (satdly) and mixing (mix) periods
    satfrq  = saturation frequency for xmtr presaturation
    satpwr  = saturation power for xmtr presaturation
    satdly  = saturation period follows D1
    slpflg  = 'y' off-resonance slp presaturation
    satshape = slpsatd (tof downfield to water)
               slpsatu (tof upfield to water)
    h2off = water off-resonance frequency (tof - satfrq)
    sspul = 'y': selects trim(x)-trim(y) at start of pulse sequence
    reverse  = 'n': gives XH proton NOESY
               'y': gives NOESY without XH protons

  phase table:

        t1 = 0 0 2 2
        t2 = 0 2 2 0 1 3 3 1 2 0 0 2 3 1 1 3 2 0 0 2 3 1 1 3 0 2 2 0 1 3 3 1
        t3 = 0 0 2 2 1 1 3 3 2 2 0 0 3 3 1 1 2 2 0 0 3 3 1 1 0 0 2 2 1 1 3 3
        t4 = 0 2 2 0 1 3 3 1 2 0 0 2 3 1 1 3 2 0 0 2 3 1 1 3 0 2 2 0 1 3 3 1
        t5 = [0 2]16
        t6 = [0 1 2 3]4
        t7 = 0 0 2 2 1 1 3 3 2 2 0 0 3 3 1 1 2 2 0 0 3 3 1 1 0 0 2 2 1 1 3 3

              DEC 12 1991   */


#include <standard.h>

pulsesequence()
{
   double          satpwr,
                   satdly,
                   satfrq,
                   mix,
                   jxh,
                   shape_pw,
                   cycle1,
                   cycle2,
                   h2off,
                   pwxlvl,
                   pwx,
                   tau;
   int             phase;
   char            sspul[MAXSTR],
                   satflg[MAXSTR],
                   slpflg[MAXSTR],
                   satshape[MAXSTR],
                   reverse[MAXSTR];

/* LOAD VARIABLES */
   satdly = getval("satdly");
   satfrq = getval("satfrq");
   satpwr = getval("satpwr");
   mix = getval("mix");
   pwxlvl = getval("pwxlvl");
   pwx = getval("pwx");
   h2off = getval("h2off");
   jxh = getval("jxh");
   tau = 1/(2.0*jxh);
   phase = (int) (getval("phase") + 0.5);
   getstr("sspul", sspul);
   getstr("slpflg",slpflg);
   getstr("satshape",satshape);
   getstr("satflg",satflg);
   getstr("reverse",reverse);

   if ((slpflg[0] == 'y') || (slpflg[2] == 'y'))
   {
    if (h2off != 0.0)
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
   if (rfwg[0] != 'y')
   {
    printf("NOTE: slpflg=y requires PPM in Channel 1. psg aborted.\n");
     abort(1);
   }
   }

/* CHECK CONDITIONS */
  
/* DETERMINE STEADY-STATE MODE */


/* STEADY-STATE PHASECYCLING

/* PHASECYCLE */

   loadtable("nosyxfil");
   
   getelem(t1,ct,v1);
   getelem(t6,ct,v7);
   assign(two,v6);

   if (reverse[0] == 'y')
    {
     getelem(t7,ct,oph);
     decr(v7);
    }
   else
    {
     getelem(t4,ct,oph);
     decr(v6);
    }


   initval(2.0*(double)(((int)(d2*getval("sw1")+0.5)%2)),v13);

     if (phase == 2)
        {incr(v1);  incr(v6);}

   add(v1,v13,v1);
   add(v6,v13,v6);
   add(oph,v13,oph);
    

/* BEGIN THE ACTUAL PULSE SEQUENCE */
   status(A);              /* relaxation status */
   rlpower(tpwr,TODEV);
   rlpower(dpwr,DODEV);
   decphase(zero);
   if (sspul[0] == 'y')
   {
      rgpulse(200*pw, zero, rof1,0.0);
      rgpulse(200*pw, one,  0.0,rof1);
   }
      hsdelay(d1);

   if (satflg[0] == 'y')
     {
        rlpower(satpwr,TODEV);
        if (slpflg[0] == 'n')
        {
        if (satfrq != tof)
         offset(satfrq,TODEV);
        rgpulse(satdly,v6,rof1,rof1);
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
          shaped_pulse(satshape,shape_pw,v6,0.0,1e-6);
         endhardloop();
         delay(rof2);
         rcvron();
         }
         else
         {
         if (satfrq != tof)
          offset(satfrq,TODEV);
         shaped_pulse(satshape,satdly,v6,rof1,rof1);
         if (satfrq != tof)
          offset(tof,TODEV);
         }
        }

        rlpower(tpwr,TODEV);
        delay(40.0e-6);
     }
   else
     delay(satdly);
   

   status(B);
     rgpulse(pw, v1, rof1, 1.0e-6);
     if (d2>0)
     delay(d2 - rof1 - 1.0e-6 - (4*pw/3.1416)); 
     else
     delay(d2);
     rgpulse(pw,t5,rof1,1.0e-6);

   status(C);
       if (satflg[2] == 'y')
        {
         rlpower(satpwr,TODEV);
         if (slpflg[2] == 'n')
         {
          if (satfrq != tof)
           offset(satfrq,TODEV);
          rgpulse(mix,v7,2.0e-6,rof1);
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
           shaped_pulse(satshape,shape_pw,v7,0.0,1e-6);
          endhardloop();
          delay(rof2);
          rcvron();
          }
          else
          {
          if (satfrq != tof)
           offset(satfrq,TODEV);
          shaped_pulse(satshape,mix,v7,2.0e-6,rof1);
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
        rgpulse(pw, t6, rof1, 0.0);
        decphase(t2);
        delay(tau);
        rcvroff();
        delay(rof1);
        decpulse(pwx,t2);
        decphase(t3);
        rgpulse(2.0*pw, t6, rof1,rof1);
        decpulse(pwx,t3);
        delay(rof1);
        rcvron();
        delay(tau-4.2e-6);
      rlpower(dpwr,DODEV);
      delay(rof2);  

    status(D);                                    /* acquisition status */
}
