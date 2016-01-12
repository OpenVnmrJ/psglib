/*  hsqc - heteronuclear single-quantum experiment using REVINEPT
           in HYPERCOMPLEX mode only 
           slp presaturation option 
           BIRD nulling option  

   pwxlvl - power level for X pulses
   pwx    - X 90 degree pulse
   jxh    - XH coupling constant
   satpwr - power for presaturation
   satdly - presaturation delay
   satfrq - presaturation frequency
   satflg - flag to turn on presaturation
   slpflg - flag to turn on slp presaturation
   satshape - slpsatd (tof downfield to water)
              slpsatu (tof upfield to water)
   h2off - water off-resonance frequency (tof - satfrq)
   null - Bird nulling delay
   sspul - 'y' turns on trim(x)-trim(y) before d1
   phase = 1,2 for hypercomplex
         ---   Krish Jan 92  */

#include <standard.h>

pulsesequence()

{
   double   pwxlvl,
            pwx,
            delta,
            bird,
            jxh,
            satpwr,
            satdly,
            satfrq,
            shape_pw,
            h2off,
            cycle1,
            cycle2,
            null,
            phase;
   int      iphase;
   char     sspul[MAXSTR],
            slpflg[MAXSTR],
            satshape[MAXSTR],
            satflg[MAXSTR];

   pwxlvl = getval("pwxlvl");
   pwx    = getval("pwx");
   jxh    = getval("jxh");
   delta  = 0.4*(1/(2*jxh));
   bird   = 1/(2*jxh);
   null   = getval("null");
   phase  = getval("phase");
   satpwr = getval("satpwr");
   satdly = getval("satdly");
   satfrq = getval("satfrq");
   h2off  = getval("h2off");

   getstr("satflg",satflg);
   getstr("slpflg",slpflg);
   getstr("satshape",satshape);
   getstr("sspul",sspul);   

   iphase = (int) (phase + 0.5);

   if ((slpflg[0] == 'y') && (h2off != 0.0))
   {
   shape_pw = 10/h2off;
   if (shape_pw < 0.0)
    shape_pw = -shape_pw;

   cycle1 = (satdly/(shape_pw + 15.4e-6));
   cycle1 = 2.0*(double)(int)(cycle1/2.0);
   initval(cycle1,v12);

   cycle2 = (null/(shape_pw + 15.4e-6));
   cycle2 = 2.0*(double)(int)(cycle2/2.0);
   initval(cycle2,v13);
   }

   if (rfwg[0] != 'y')
   {
    if ((slpflg[0] == 'y') || (slpflg[1] == 'y'))
    {
     printf("NOTE: slpflg=y requires PPM in channel 1. PSG aborted.\n");
     abort(1);
    }
   }

   loadtable("slphsqc");          /* t1 = 1 1 3 3
                                  t2 = 0 2
                                  t3 = 0 0 0 0 2 2 2 2
                                  t4 = 0 0 0 0 0 0 0 0 2 2 2 2 2 2 2 2
                                  t5 = 0 2 2 0 2 0 0 2 2 0 0 2 0 2 2 0
                               */

   getelem(t2,ct,v2);
   getelem(t5,ct,oph);

   initval(2.0*(double)(((int)(d2*getval("sw1")+0.5)%2)),v14);

   if (iphase == 2)
     incr(v2);

   add(v2,v14,v2);
   add(oph,v14,oph);

   status(A);
     rlpower(pwxlvl,DODEV);
     if (sspul[0] == 'y')
      {
       rgpulse(200*pw,zero,rof1,0.0);
       rgpulse(200*pw,one,0.0,rof1);
      }
      hsdelay(d1);

     if (satflg[0] == 'y')
      {
       rlpower(satpwr,TODEV);
       if (slpflg[0] == 'y')
       {
        if (h2off != 0.0)
        {
        rcvroff();
        delay(rof1);
        starthardloop(v12);
         shaped_pulse(satshape,shape_pw,zero,0.0,1e-6);
        endhardloop();
        delay(rof2);
        rcvron();
        }
        else
        {
        if (satfrq != tof)
         offset(satfrq,TODEV);
        shaped_pulse(satshape,satdly,zero,4.0e-5, 0.2e-6);
        if (satfrq != tof)
         offset(tof,TODEV);
        }
       }
       else
       {
        if (satfrq != tof)
         offset(satfrq,TODEV);
        rgpulse(satdly,zero,4.0e-5, 0.2e-6);
        if (satfrq != tof)
         offset(tof,TODEV);
       }
       rlpower(tpwr,TODEV);
       delay(1.0e-5);
      }

    status(B);
     if (null != 0.0)
      {
       rcvroff();
       rgpulse(pw,zero,1.0e-5,0.0);
       delay(bird);
       decrgpulse(pwx,one,rof1,0.0);
       simpulse(2.0*pw,2.0*pwx,zero,zero,1.0e-6,0.0);
       decphase(one);
       decrgpulse(pwx,one,1.0e-6,0.0);
       delay(bird);
       rgpulse(pw,two,rof1,rof2);
       rcvron();
       if (satflg[1] == 'y')
        {
         rlpower(satpwr,TODEV);
         if (slpflg[1] == 'y')
         {
          if (h2off != 0.0)
          {
          rcvroff();
          delay(rof1);
          starthardloop(v13);
           shaped_pulse(satshape,shape_pw,zero,0.0,1e-6);
          endhardloop();
          delay(rof2);
          rcvron();
          }
          else
          {
          if (satfrq != tof)
           offset(satfrq,TODEV);
          shaped_pulse(satshape,null,zero,4.0e-5,1.0e-5);
          if (satfrq !=tof)
           offset(tof,TODEV);
          }
         }
         else
         {
          if (satfrq != tof)
           offset(satfrq,TODEV);
          rgpulse(null,zero,4.0e-5,1.0e-5);
          if (satfrq !=tof)
           offset(tof,TODEV);
         }
         rlpower(tpwr,TODEV);
         delay(4.0e-5);
        }
       else
        hsdelay(null);
      }
     txphase(zero);
     rcvroff();
     rgpulse(pw,zero,rof1,0.0);
     decphase(zero);
     txphase(zero);
     delay(delta - pwx);
     simpulse(2*pw,2*pwx,zero,zero,0.0,0.0);
     txphase(t1);
     decphase(v2);
     delay(delta - pwx);
     simpulse(pw,pwx,t1,v2,0.0,0.0);
     txphase(zero);
     decphase(t4);
     delay(d2/2);
     rgpulse(2*pw,zero,0.0,0.0);
     txphase(t3);
     decphase(t4);
     delay(d2/2);
     simpulse(pw,pwx,t3,t4,0.0,0.0);
     txphase(zero);
     decphase(zero);
     delay(delta - pwx - (2*pwx/3.1416));
     simpulse(2*pw,2*pwx,zero,zero,0.0,0.0);
     rlpower(dpwr,DODEV);
     delay(rof2);
     rcvron();
     delay(delta - pwx - 4.2e-6 - rof2);
   status(C);
}
