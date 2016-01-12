/* bashDnosy -  BAnd-Selected Homonuclear Decoupled NOESY 
            States-Hypercomplex only

	Note: very similar to bashdtoxy.  See bashdtoxy for more details.
	Krish	Aug, 1995	*/


#include <standard.h>

static int	ph1[2] = {0,2},
		ph2[1] = {0},
		ph6[2] = {0,2},
		ph3[8] = {0,2,1,3,2,0,3,1},
		ph4[8] = {0,0,1,1,2,2,3,3};

pulsesequence()
{
   double	   hsgpwr,
		   mix,
		   gzlvl1,
		   gt1,
		   gzlvl2,
		   gt2,
		   gstab,
		   selpwr,
		   alfa1,
		   t1dly,
		   selpw;
   int             iphase;
   char            sspul[MAXSTR],
		   pwshape[MAXSTR],	
		   fadflg[MAXSTR],
		   homodec[MAXSTR],
                   satflg[MAXSTR];

   hsgpwr = getval("hsgpwr");
   gzlvl1 = getval("gzlvl1");
   gt1 = getval("gt1");
   gzlvl2 = getval("gzlvl2");
   gt2 = getval("gt2");
   gstab =getval("gstab");
   selpwr = getval("selpwr");
   selpw = getval("selpw");
   mix = getval("mix");
   getstr("pwshape",pwshape);
   getstr("fadflg",fadflg);
   getstr("homodec",homodec);
   iphase = (int) (getval("phase") + 0.5);
   getstr("sspul", sspul);
   getstr("satflg",satflg);

   alfa1 = 2.0e-6+(4*pw/PI);
   if (homodec[0] == 'y')
	alfa1 = alfa1 + 2.0e-6 + 2*pw;
   t1dly = d2-alfa1;
   if (t1dly > 0.0)
	t1dly = t1dly;
   else
	t1dly = 0.0;

   settable(t1,2,ph1);
   settable(t4,8,ph4);
   settable(t2,1,ph2);
   settable(t6,2,ph6);
   settable(t3,8,ph3);

   getelem(t1,ct,v1);
   getelem(t4,ct,oph);
   assign(zero,v6);

   if (iphase == 2) 
      {incr(v1);  incr(v6); }

   assign(v1,v8);
   add(v8,two,v7);

   initval(2.0*(double)(((int)(d2*getval("sw1")+0.5)%2)),v14);
   if (fadflg[0] == 'y')
   {
   add(v1, v14, v1);
   add(oph,v14,oph);
   add(v6,v14,v6);  
   }

/* BEGIN THE ACTUAL PULSE SEQUENCE */
   status(A);
   obspower(tpwr);
   delay(5.0e-5);

   if (sspul[0] == 'y')
   {
	zgradpulse(hsgpwr,2.0e-3);
	rgpulse(pw,zero,rof1,rof1);
	zgradpulse(hsgpwr,2.0e-3);
   }

   delay(d1);

   if (satflg[0] == 'y') 
     {
        obspower(satpwr);
	if (satfrq != tof)
	obsoffset(satfrq);
	rgpulse(satdly,v6,rof1,rof1);
	if (satfrq != tof)
	obsoffset(tof);
        obspower(tpwr);
     }

   status(B);
      rgpulse(pw, v1, rof1, 1.0e-6);

      if (homodec[0] == 'y')
       {
                delay(t1dly/2);
                rgpulse(2*pw,t6,1.0e-6,1.0e-6);
       }
      else
                delay(t1dly);
	delay(gstab);
        zgradpulse(gzlvl1,gt1);
        delay(gstab);
        obspower(selpwr);
        shaped_pulse(pwshape,selpw,v8,rof1,rof1); 
        obspower(tpwr);
        delay(gstab); 
        zgradpulse(gzlvl1,gt1);
        delay(gstab);
      if (homodec[0] == 'y')
                delay(t1dly/2);
        delay(gstab); 
        zgradpulse(gzlvl2,gt2);
        delay(gstab); 
        obspower(selpwr); 
        shaped_pulse(pwshape,selpw,v7,rof1,rof1);  
        obspower(tpwr); 
        delay(gstab); 
        zgradpulse(gzlvl2,gt2);
        delay(gstab);

   status(C);

      rgpulse(pw,t2,1.0e-6,1.0e-6); 
   if (satflg[2] == 'y')
     {
        obspower(satpwr);
	if (satfrq != tof)
	obsoffset(satfrq);
	rgpulse(mix,zero,rof1,rof1);
	if (satfrq != tof)
	obsoffset(tof);
        obspower(tpwr);
     }
   else
      delay(mix);

      rgpulse(pw,t3,1.0e-6,rof2);

   status(D);
}
