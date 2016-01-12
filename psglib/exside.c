/*  exside - Excitation Sculptured Indirect Detection Experiment
		J-scaling option

	processing : wft2d(1,0,1,0,0,-1,0,1)
		rp = rp(s2pul) - 90
		rp1 = 0

krish krishnamurthy Aug. 95

[REF: V.V. Krishnamurthy, JMR, A121, 33-41 (1996)]
*/

#include <standard.h>

static int	ph4[8]  =  {0,0,0,0,2,2,2,2},
		ph2[2]  =  {0,2},
		ph1[4]  =  {0,0,2,2},
		ph5[8]  =  {0,2,2,0,2,0,0,2};

pulsesequence()

{
   double   hsgt,
	    gzlvl1,
	    gt1,
	    gzlvl2,
	    gt2,
	    gzlvl3,
	    gt3,
	    gzlvl4,
	    gt4,
	    t1dly,
	    jdly,
	    gstab,
	    bird,
	    hsgpwr,
	    selpwr,
	    selpw;
   int      iphase,
	    ijscale,
	    icosel;
   char     sspul[MAXSTR],
	    nullflg[MAXSTR],
	    pwshape[MAXSTR];

   hsgt = getval("hsgt");
   gzlvl1 = getval("gzlvl1");
   gt1 = getval("gt1");
   gzlvl2 = getval("gzlvl2");
   gt2 = getval("gt2");
   gzlvl3 = getval("gzlvl3");
   gt3 = getval("gt3");
   gzlvl4 = getval("gzlvl4");
   gt4 = getval("gt4");
   gstab = getval("gstab");
   hsgpwr = getval("hsgpwr");
   selpwr = getval("selpwr");
   selpw = getval("selpw");
   tau   = 1/(4*(getval("jnxh")));
   bird = 1/(4*(getval("j1xh")));
   getstr("sspul",sspul);
   getstr("nullflg",nullflg);
   getstr("pwshape",pwshape);
   ijscale = (int)(getval("jscale") + 0.5);

   t1dly = d2 - 2*pw - 2.0e-6;
   if (t1dly < 0.0) t1dly = 0.0;
   jdly = (ijscale*d2);
   if (jdly < 0.0) jdly = 0.0;

   iphase = (int) (getval("phase") + 0.5);

   settable(t1,4,ph1);
   settable(t2,2,ph2);
   settable(t4,8,ph4);
   settable(t5,8,ph5);

   getelem(t2,ct,v2);
   getelem(t5,ct,oph);

   initval(2.0*(double)(((int)(d2*getval("sw1")+0.5)%2)),v14);

   if (iphase == 2)
     icosel = -1;
   else
     icosel = 1;

   add(v2,v14,v2);
   add(oph,v14,oph);

   status(A);
     decpower(pwxlvl);
     obspower(tpwr);
     delay(5.0e-5);

     if (sspul[0] == 'y')
      {
	zgradpulse(hsgpwr,0.002);
	rgpulse(pw,zero,rof1,rof1);
	zgradpulse(hsgpwr,0.002);
      }

      delay(d1);
       if (nullflg[0] == 'y')
         {
          rgpulse(0.5*pw,zero,rof1,rof1);
          delay(2*bird);
          simpulse(2.0*pw,2.0*pwx,zero,zero,rof1,rof1);
          delay(2*bird);
          rgpulse(1.5*pw,zero,rof1,rof1);
          zgradpulse(hsgpwr,0.01);
          delay(0.05);
         }

    status(B);
     rcvroff();
     rgpulse(pw,t1,rof1,rof1);
     delay(jdly/4);
     delay(tau/2 - gt1 - gstab);
     zgradpulse(gzlvl1,gt1);
     delay(gstab+WFG_STOP_DELAY);
     obspower(selpwr);
     simshaped_pulse(pwshape,"hard",selpw,2*pwx,zero,zero,rof1,rof1);
     obspower(tpwr);
     delay(gstab+WFG_START_DELAY);
     zgradpulse(gzlvl1,gt1);
     delay(tau/2 - gt1 - gstab);
     delay(jdly/2);
     delay(tau/2 - gt2 - gstab);
     zgradpulse(gzlvl2,gt2);
     delay(gstab+WFG_STOP_DELAY);
     obspower(selpwr);
     simshaped_pulse(pwshape,"hard",selpw,2*pwx,two,zero,rof1,rof1);
     obspower(tpwr);
     delay(gstab+WFG_START_DELAY);
     zgradpulse(gzlvl2,gt2);
     delay(tau/2 - gt2 - gstab);
     delay(jdly/4);
     rgpulse(pw,one,rof1,rof1);
     zgradpulse(hsgpwr,hsgt);
     delay(1.0e-4);
     decrgpulse(pwx,v2,rof1,1.0e-6);

     delay(t1dly/2);
     rgpulse(2*pw,zero,1.0e-6,1.0e-6);
     delay(t1dly/2);

     delay(gt3+gstab+gt3+GRADIENT_DELAY);
     decrgpulse(2*pwx,v2,1.0e-6,1.0e-6);
     delay(gstab);
     zgradpulse(gzlvl3,gt3);
     delay(gt3);
     simpulse(pw,pwx,zero,t4,1.0e-6,1.0e-6);
     delay(gstab);
     zgradpulse(icosel*gzlvl4,gt4);
     delay(tau/2 - gt1 - gstab - gt4 - gstab - GRADIENT_DELAY);
     zgradpulse(gzlvl1,gt1);
     delay(gstab+WFG_STOP_DELAY);
     obspower(selpwr);
     simshaped_pulse(pwshape,"hard",selpw,2*pwx,zero,zero,rof1,rof1);
     obspower(tpwr);
     delay(gstab+WFG_START_DELAY);
     zgradpulse(gzlvl1,gt1);
     delay(tau - gt1 - gt2 - 2*gstab);
     zgradpulse(gzlvl2,gt2);
     delay(gstab+WFG_STOP_DELAY);
     obspower(selpwr);
     simshaped_pulse(pwshape,"hard",selpw,2*pwx,two,zero,rof1,rof2);
     obspower(tpwr);
     delay(gstab+WFG_START_DELAY);
     zgradpulse(gzlvl2,gt2);
     delay(tau/2 - gt2 - gstab - POWER_DELAY);
     decpower(dpwr);
     rcvron();

   status(C);
}
