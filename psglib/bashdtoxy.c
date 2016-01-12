/* bashDtoxy -  BAnd-Selected Homonuclear Decoupled TOCSY 
            Rance type Sensitivity enhanced

	Band selection is based on double PFG spin echo


     homodec = 'y'
	wft2d(1,0,0,-1,-1,0,0,-1,0,-1,-1,0,0,-1,1,0)
		small non-zero rp1 and lp1 is expected

     homodec = 'n'
        wft2d(1,0,0,1,-1,0,0,1,0,1,-1,0,0,1,1,0)
		rp1=180 lp1=0

	Krish Krishnamurthy	Aug, 1995
	[REF: V.V. krishnamurthy, Magn. Reson. Chem., in press (as of Aug. 1996)]
	*/


#include <standard.h>
extern int dps_flag;

mleva()
{
 double slpw;
 slpw = getval("slpw");
 txphase(v2); delay(slpw); 
 xmtroff(); delay(slpw); xmtron();
 txphase(v3); delay(2*slpw); 
 xmtroff(); delay(slpw); xmtron();
 txphase(v2); delay(slpw);
}

mlevb()
{
 double slpw;
 slpw = getval("slpw");
 txphase(v4); delay(slpw); 
 xmtroff(); delay(slpw); xmtron();
 txphase(v5); delay(2*slpw); 
 xmtroff(); delay(slpw); xmtron();
 txphase(v4); delay(slpw);
}

dipsi(phs1,phs2)
codeint phs1,phs2;
{
	double slpw5;
	slpw5 = getval("slpw")/18.0;

	txphase(phs1); delay(64*slpw5);
	txphase(phs2); delay(82*slpw5);
	txphase(phs1); delay(58*slpw5);
	txphase(phs2); delay(57*slpw5);
	txphase(phs1); delay(6*slpw5);
	txphase(phs2); delay(49*slpw5);
	txphase(phs1); delay(75*slpw5);
	txphase(phs2); delay(53*slpw5);
	txphase(phs1); delay(74*slpw5);
}

static int	ph1[8] = {0,2,0,2,1,3,1,3},
		ph2[8] = {1,1,3,3,2,2,0,0},
		ph3[8] = {0,0,0,0,1,1,1,1},
		ph5[8] = {0,0,0,0,1,1,1,1},
		ph6[8] = {0,2,0,2,1,3,1,3},
		ph4[8] = {0,2,0,2,1,3,1,3};

pulsesequence()
{
   double	   hsgpwr,
		   slpwr,
		   slpw,
		   slpw5,
		   slmix,
		   cycles,
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
		   dipsiflg[MAXSTR],
		   homodec[MAXSTR],
		   fadflg[MAXSTR],
                   satflg[MAXSTR];

   hsgpwr = getval("hsgpwr");
   slpwr = getval("slpwr");
   slmix = getval("slmix");
   slpw = getval("slpw");
   slpw5 = slpw/18.0;
   gzlvl1 = getval("gzlvl1");
   gt1 = getval("gt1");
   gzlvl2 = getval("gzlvl2");
   gt2 = getval("gt2");
   gstab =getval("gstab");
   selpwr = getval("selpwr");
   selpw = getval("selpw");
   getstr("pwshape",pwshape);
   getstr("dipsiflg",dipsiflg);
   getstr("homodec",homodec);
   getstr("fadflg",fadflg);
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

   if (dipsiflg[0] == 'y')
     cycles = slmix/(2072*slpw5);
   else
     cycles = (slmix)/(96.67*slpw);
   cycles = (double)(int)(cycles);
   initval(cycles,v9);

   settable(t1,8,ph1);
   settable(t3,8,ph3);
   settable(t5,8,ph5);
   settable(t4,8,ph4);
   settable(t6,8,ph6);
   settable(t2,8,ph2);

   getelem(t1,ct,v1);
   getelem(t4,ct,oph);
   getelem(t3,ct,v10);
   getelem(t2,ct,v2);
   assign(zero,v6);

   add(v2,one,v3);
   add(v3,one,v4);
   add(v4,one,v5);

   if (iphase == 2) 
      {incr(v1);  incr(v6); }
   if (iphase == 3)
      add(v10,two,v10);
   if (iphase == 4)
       {incr(v1); incr(v6); add(v10,two,v10);}

   assign(v1,v8);
   add(v8,two,v7);

   initval(2.0*(double)(((int)(d2*getval("sw1")+0.5)%2)),v14);
  if (fadflg[0] == 'y')
  {
   add(v1, v14, v1);
   add(v6,v14,v6);
   add(oph,v14,oph);
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
	rgpulse(satdly,v6,rof1,rof2);
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
	rcvroff();
	rgpulse(pw,v10,1.0e-6,0.0);
      obspower(slpwr);
      if (cycles > 0.0)
      {
        if (dps_flag)
	 rgpulse(slmix,zero,rof1,rof1);
	else
	{
	xmtron();
	if (dipsiflg[0] == 'y')
	 {
		starthardloop(v9);
                dipsi(v2,v4);
                dipsi(v4,v2);
                dipsi(v4,v2);
                dipsi(v2,v4);
		endhardloop();
	 }
	else
	 {
          starthardloop(v9);
            mleva(); mleva(); mlevb(); mlevb();
            mlevb(); mleva(); mleva(); mlevb();
            mlevb(); mlevb(); mleva(); mleva();
            mleva(); mlevb(); mlevb(); mleva();
            txphase(v3); delay(0.67*slpw);
          endhardloop();
	 }
        xmtroff();
	}
       }
      obspower(tpwr);
      rgpulse(pw,t5,0.0,rof2);
      rcvron();

   status(D);
}
