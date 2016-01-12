/* trnsySS3d - A sequence to do either	TOCSY-NOESY-3D
					ROESY-NOESY-3D
					TOCSY-ROESY-3D
		with SS read pulse.
		(For tocsy-roesy: Z-filter before SS read pulse)

		Gradient only during the NOESY mixing time/z-filter delay  

  Parameters:
	General:
		nosyflg, toxyflg and rosyflg - flags to turn on/off
		SSpwr - power for SS pulse
		pwSS - SS pulse width
		SSshape - SS shape
		pw0 - "stearing" pulse (same phase as pwSS) - can be negative
		pw1 - "stearing" pulse (orthogonal to pwSS) - can be negative
		d2corr, d3corr - enterable evolution delay corrections
				    Can be negative.
		gt - gradient time
		gstab - time for system recovery after gradient
		gzlvl - gradient level
		zfilt - z-filter time (zfilt must be > gstab+gt)

	TOCSY:
		slpwrT - spinlock power
		slpwT  - pulse width for mlev17
		wdwfctr - window factor in mlev17 (window = wfwfctr*slpwT)
		mixT - tocsy mixing time
		trim - trim pulse

	ROESY:
		slpwrR - spinlock power (typically same as tpwr)
		slpwR - pulse width for time-shared spinlock (typically
					30 degree pulse as in roesy)
		ratio - (similar to that in roesy)
		xtradly - NOT ENTERABLE (correction for lp=0
					in roesy evolution axis)
		mixR - roesy mixing time

	NOESY:
		mixN - noesy mixing time
 
   Processing:

	2D processing:

	array: phase,phase2
		wft2d('ni',#,1,0,0,0,0,0,0,0,0,0,0,0,-1,0,0,0)
		wft2d('ni2',#,1,0,0,0,0,0,0,0,0,0,-1,0,0,0,0,0)

	array: phase2,phase
		wft2d('ni',#,1,0,0,0,0,0,0,0,0,0,-1,0,0,0,0,0)
		wft2d('ni2',#,1,0,0,0,0,0,0,0,0,0,0,0,-1,0,0,0)

	3D processing with N-type selection!!



   NOTE:   For ROESY-NOESY, roesy acquisition time (i.e., noesy evolution
		time) must be corrected to make lp2=0.  Use a macro similar
		to calfa to back calculate the correction (d3corr in this
		case) and enter it in the parameter.   Typically 5-10 us
		correction may be needed.  Roesy evolution doesn't seem
              	to require correction!! [See xtradly in the sequence]
           For TOCSY-ROESY, Roesy evolution doesn't seem to require 
		correction!! [See xtradly in the sequence].
		alfa must be corrected as usual using calfa.
	   For TOCSY-NOESY,  there seems no need for corrections.


 
Krish Krishnamurthy (Sept 27, 1993)

      */


#include <standard.h>
#define PI 3.1416

mleva()
{
  double window,slpwT;
  slpwT=getval("slpwT");
  window=getval("wdwfctr")*slpwT;
  rgpulse(slpwT,v5,0.0,window);
  rgpulse(2*slpwT,v6,0.0,window);
  rgpulse(slpwT,v5,0.0,0.0);
}

mlevb()
{
  double window,slpwT;
  slpwT = getval("slpwT");
  window=getval("wdwfctr")*slpwT;
  rgpulse(slpwT,v7,0.0,window);
  rgpulse(2*slpwT,v8,0.0,window);
  rgpulse(slpwT,v7,0.0,0.0);
}

static int 	ph1[2] = {0,2},
		ph2[1] = {0},
		ph3[8] = {0,2,2,0,1,3,3,1},
		ph4[16] = {1,1,1,1,1,1,1,1,3,3,3,3,3,3,3,3},
		ph5[16] = {0,0,0,0,0,0,0,0,2,2,2,2,2,2,2,2},
		ph6[8] = {0,0,2,2,1,1,3,3},
                ph7[8] = {0,2,2,0,1,3,3,1},
		ph8[8] = {0,2,2,0,1,3,3,1};

/* ph1 - first pulse
   ph2 - flip back (i.e., 2nd pulse in noesy)
   ph3 - SS pulse
   ph4 - roesy spinlock
   ph5 - tocsy mlev17 spinlock
   ph6 - receiver
   ph7 - first stearing pulse
   ph8 - second stearing pulse            */


pulsesequence()
{
   double 	phase,
		phase2,
		pwSS,
		SSpwr,
		pw0,
		pw1,
		mixT,
                mixR,
                mixN,
		slpwrT,
		slpwT,
                trim,
                wdwfctr,
                window,
                slpwrR,
		slpwR,
		ratio,
		xtradly,
		d2corr,
		d3corr,
		cyclesT,
		cyclesR,
                ss,
		gstab,
		gt,
		gzlvl,
		zfilt;
   char 	sspul[MAXSTR], 
		SSshape[MAXSTR],
		toxyflg[MAXSTR],
		rosyflg[MAXSTR],
		nosyflg[MAXSTR];
   int		iphase,
		iphase2;

   ss = getval("ss");
   mixT=getval("mixT");
   mixR=getval("mixR");
   mixN=getval("mixN");
   slpwrT=getval("slpwrT");
   slpwT=getval("slpwT");
   slpwrR=getval("slpwrR");
   slpwR=getval("slpwR");
   wdwfctr=getval("wdwfctr");
   window=(wdwfctr*slpwT);
   trim = getval("trim");
   ratio = getval("ratio");
   d2corr = getval("d2corr");
   d3corr = getval("d3corr");
   gstab=getval("gstab");
   gt=getval("gt");
   gzlvl=getval("gzlvl");
   zfilt=getval("zfilt");
   cyclesR = mixR/(16.0*(ratio+1.0)*slpwR);
   cyclesR = 2.0*(double)(int)(cyclesR/2.0);
   initval(cyclesR,v10);
   cyclesT = (mixT - trim)/(64.66*slpwT+32*window);
   cyclesT = 2.0*(double)(int)(cyclesT/2.0);
   initval(cyclesT,v9);
   pwSS=getval("pwSS");
   pw0=getval("pw0");
   pw1=getval("pw1");
   SSpwr = getval("SSpwr");
   getstr("sspul",sspul);
   getstr("SSshape",SSshape);
   getstr("toxyflg",toxyflg);
   getstr("rosyflg",rosyflg);
   getstr("nosyflg",nosyflg);
   phase = getval("phase");
   phase2 = getval("phase2");
   iphase = (int)(phase+0.5);
   iphase2 = (int)(phase2+0.5);
   xtradly = (ratio-1)*8e-6;
   xtradly = xtradly - POWER_DELAY;
   mixN = (mixN - zfilt);
   if (mixN < 0.0)
    mixN = 0.0;

   if (toxyflg[0] == 'y')
   {
     if ((rosyflg[0] == 'y') && (nosyflg[0] == 'y'))
      {
        printf("ALL THREE FLAGS ARE y !!! - ACQUISITION ABORTED.\n");
        abort(1);
      }
   }


   initval(7.0,v2);
   stepsize(45.0,TODEV);

   if (ss < 0)
	ss = (-1)*ss;
   else
	{  if ((ss>0) && (ix == 1))
		ss=ss;
	   else
		ss=0;
	}
   initval(ss,ssctr);
   initval(ss,ssval);

   ifzero(ssctr);
	assign(ct,v11);
   elsenz(ssctr);
	sub(ssval,ssctr,v11);
   endif(ssctr);

   settable(t1,2,ph1);
   settable(t2,1,ph2);
   settable(t3,8,ph3);
   settable(t4,16,ph4);
   settable(t5,16,ph5);
   settable(t6,8,ph6);
   settable(t7,8,ph7);
   settable(t8,8,ph8);

    getelem(t1,v11,v1);
    getelem(t3,v11,v12);
    getelem(t4,v11,v4);
    getelem(t5,v11,v5);
    getelem(t6,v11,oph);

    if (iphase == 2)
     incr(v1);
    if (iphase2 == 2)
    { incr(v1);
      if (toxyflg[0]=='y')
       incr(v5);
      else
       incr(v4);
     }

    initval(2.0*(double)(((int)(d2*getval("sw1")+0.5)%2)),v13);
    initval(2.0*(double)(((int)(d3*getval("sw2")+0.5)%2)),v14);
    add(v1,v14,v1);
    add(oph,v14,oph);
    add(v1,v13,v1);
    add(oph,v13,oph);
    if (toxyflg[0] == 'y')
     add(v5,v14,v5);
    else
     add(v4,v14,v4);

    add(v5,one,v6);
    add(v6,one,v7);
    add(v7,one,v8);

    if (pw0 < 0.0)
    {  tsadd(t7,2,4);
       pw0=(-1.0)*pw0;}

    if (pw1 < 0.0)
    {  add(t8,2,4);
       pw1=(-1.0)*pw1;}

   status(A);
     rlpower(tpwr,TODEV);
     if (sspul[A] == 'y')
      { rlpower(tpwr-12,TODEV);
        rgpulse(200*pw,zero,rof1,0.0);
        rgpulse(200*pw,one,0.0,rof1);
        rlpower(tpwr,TODEV); }
     xmtrphase(v2);
     txphase(v1);
     hsdelay(d1);
     rcvroff();
     rgpulse(pw, v1, rof1, 1.0e-6);
     xmtrphase(zero);

   status(B);
      if (toxyflg[0] == 'y')
      {
	rlpower(slpwrT,TODEV);
	if (d2 > 0.0)
	  delay(d2 - POWER_DELAY - 1.0e-6 - (2*pw/PI) + d2corr);
        else
          delay(d2 + d2corr);
        if (cyclesT > 1.0)
         {
           rgpulse(trim,v6,0.0,0.0);
           starthardloop(v9);
		mleva(); mleva(); mlevb(); mlevb();
		mlevb(); mleva(); mleva(); mlevb();
		mlevb(); mlevb(); mleva(); mleva();
		mleva(); mlevb(); mlevb(); mleva();
		rgpulse(0.66*slpwT,v6,0.0,0.0);
	   endhardloop();
         }
        rlpower(tpwr,TODEV);
      }

    status(C);
      if (rosyflg[0] == 'y')
      {
       if (toxyflg[0] == 'y')
       {
        if (d3>0.0)
         delay(d3 - (2*pw/PI) - POWER_DELAY - rof1 + d3corr);
        else
         delay(d3 + d3corr);
       }
       else
       {
        if (d2 > 0.0)
         delay(d2 - (4*pw/PI) - 1.0e-6 - rof1 + d2corr);
        else
         delay(d2 + d2corr);
       }
       delay(xtradly);
       rgpulse(pw,v4,rof1,0.0);
       rlpower(slpwrR,TODEV);
       delay(ratio*slpwR/2);
       if (cyclesR > 1.0)
       {
       starthardloop(v10);
	rgpulse(slpwR,v4,0.0,ratio*slpwR);
	rgpulse(slpwR,v4,0.0,ratio*slpwR);
	rgpulse(slpwR,v4,0.0,ratio*slpwR);
	rgpulse(slpwR,v4,0.0,ratio*slpwR);
	rgpulse(slpwR,v4,0.0,ratio*slpwR);
	rgpulse(slpwR,v4,0.0,ratio*slpwR);
	rgpulse(slpwR,v4,0.0,ratio*slpwR);
	rgpulse(slpwR,v4,0.0,ratio*slpwR);
	rgpulse(slpwR,v4,0.0,ratio*slpwR);
	rgpulse(slpwR,v4,0.0,ratio*slpwR);
	rgpulse(slpwR,v4,0.0,ratio*slpwR);
	rgpulse(slpwR,v4,0.0,ratio*slpwR);
	rgpulse(slpwR,v4,0.0,ratio*slpwR);
	rgpulse(slpwR,v4,0.0,ratio*slpwR);
	rgpulse(slpwR,v4,0.0,ratio*slpwR);
	rgpulse(slpwR,v4,0.0,ratio*slpwR);
       endhardloop();
       }
       rlpower(tpwr,TODEV);
       delay(ratio*slpwR/2);
       rgpulse(pw,v4,0.0,0.0);
     }

   status(D);
      if (nosyflg[0] == 'y')
      {
        if (toxyflg[0] == 'y')
        { if (d3 > 0)
          delay(d3 -  POWER_DELAY - (2*pw/PI) - rof1 + d3corr);
          else
          delay(d3 + d3corr);
        }
        else
        {
          delay(xtradly);
          if (d3 > 0)
          delay(d3 - (4*pw/PI) - rof1 + d3corr);
          else
          delay(d3 + d3corr);
        }
      }

   status(E);
      rgpulse(pw, t2, rof1, 1.0e-6);
      if (nosyflg[0] == 'y')
       delay(mixN - zfilt);
      delay(zfilt-gt-gstab);
      zgradpulse(gzlvl,gt);
      delay(gstab);
      rlpower(SSpwr,TODEV);
      shaped_pulse(SSshape,pwSS,v12,rof1,0.0);
      rlpower(SSpwr-12,TODEV);
      rgpulse(pw0,t7,rof1,0.0);
      rgpulse(pw1,t8,0.0,rof1);
      rcvron();

   status(F);
}
