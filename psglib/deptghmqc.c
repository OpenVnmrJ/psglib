/* depthmqc.c  RC version 1.0 JMR 85,400-405(1989) */

/* modifief BKJ 8 July 96 */

#include <standard.h>

static int	phs1[8] = {0,0,0,0,2,2,2,2},
		phs4[2] = {0,2},
		phs3[8] = {1,1,1,1,3,3,3,3},
		phs2[4] = {0,0,2,2},
		phs5[8] = {0,2,2,0,2,0,0,2};

pulsesequence()
{
	double	pwx, 
		pwxlvl, 
		mult, 
		hsgpwr,
		hsgt,
		t1dly,
		gzlvl1,
		gt1,
		gstab,
		gzlvl3,
		gt3,
		gtau,
		tauc,
		alfa1,
		bird;
	int	iphase,
		icosel;
   	char	sspul[MAXSTR],
		nullflg[MAXSTR];
		

pwx = getval("pwx");
pwxlvl = getval("pwxlvl");
mult = getval("mult");
getstr("sspul",sspul);
gzlvl1 = getval("gzlvl1");
gt1 = getval("gt1");
gstab=getval("gstab");
gzlvl3 = getval("gzlvl3");
gt3 = getval("gt3");
alfa1 = getval("alfa1");
getstr("nullflg",nullflg);
hsgpwr = getval("hsgpwr");
hsgt = getval("hsgt");
iphase = (int)(getval("phase") + 0.5);
bird = 1/(2 * (getval("j1xh")));
alfa1 = alfa1 + 4.0e-6 + 2*pw;
gtau = 2*gstab + GRADIENT_DELAY;
tauc = gtau + gt3 - gstab;
t1dly = d2 - alfa1;
if (t1dly < 0.0)
 t1dly = 0.0;

settable(t1,8,phs1);
settable(t2,4,phs2);
settable(t3,8,phs3);
settable(t4,2,phs4);
settable(t5,8,phs5);

getelem(t5, ct, oph);
getelem(t2, ct, v2);   

assign(zero,v4);

if (iphase == 2)
	icosel = 1;
else
	icosel = -1;

initval(2.0*(double)(((int)(d2*getval("sw1")+0.5)%2)),v14);

add(v2,v14,v2);
add(v4,v14,v4);
add(oph,v14,oph);


status(A);
decpower(pwxlvl);
obspower(tpwr);
rcvroff();
if (sspul[0] == 'y')
	{
		zgradpulse(hsgpwr,0.01);
		rgpulse(pw,zero,rof1,rof1);
                zgradpulse(hsgpwr,0.01); 
	}
delay(d1); 

status(B);
if (nullflg[0] == 'y')
	{
		rgpulse(0.5*pw, zero, rof1, 2.0e-6);
		delay(bird);
		simpulse(2*pw, 2*pwx, zero,zero, 2.0e-6, 2.0e-6);
		delay(bird);
		rgpulse(1.5*pw, two, 2.0e-6, rof1);
		zgradpulse(hsgpwr,hsgt);
		delay(0.01);
	}

rgpulse(pw, t1, rof1, 2.0e-6);
delay(bird);
simpulse(2*pw, pwx, zero, v2, 2.0e-6, 2.0e-6);
delay(bird);
simpulse(mult*pw, 2*pwx, t3, v4, 2.0e-6, 2.0e-6);
delay(bird);

delay(gt1+gtau);
decrgpulse(2*pwx,v4,2.0e-6,2.0e-6);
delay(gstab);
zgradpulse(gzlvl1,gt1);
delay(gstab);

delay(t1dly/2.0);
rgpulse(2*pw, t1, 2.0e-6,2.0e-6);
delay(t1dly/2.0);

delay(gstab);
zgradpulse(gzlvl1,gt1);
delay(gstab);
decrgpulse(2*pwx,zero,2.0e-6,2.0e-6);
delay(gt1+gtau - (2*pwx/PI));
decrgpulse(pwx, t4, 2.0e-6,2.0e-6);
delay(gstab);
zgradpulse(icosel*gzlvl3,gt3);
decpower(dpwr);
delay(bird - tauc);
delay(rof2 - POWER_DELAY);
rcvron();
status(C);
}

