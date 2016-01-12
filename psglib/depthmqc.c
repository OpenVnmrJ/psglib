/* depthmqc.c  RC version 1.0 JMR 85,400-405(1989) */

#include <standard.h>

static double	d2_init = 0.0;

pulsesequence()
{
	double	pwx, pwxlvl, j, mult, bird,
			null, sw1, sstrim;
	int		iphase, t1_counter;

pwx = getval("pwx");
pwxlvl = getval("pwxlvl");
j = getval("j");
mult = getval("mult");
null = getval("null");
sstrim = getval("sstrim");
sw1 = getval("sw1");
iphase = (int)(getval("phase") + 0.5);

if (j != 0.0)
	bird = 1.0/(2.0*j);
else
	bird = 3.1e-3;

loadtable("depthmqc");
/*		t1 = 0 0 0 0 2 2 2 2
		t2 = 0 2 0 2
		t3 = 1 1 1 1 3 3 3 3
		t4 = 0 0 2 2
		t5 = 0 2 2 0 2 0 0 2
		t8 = 1 1 1 1 3 3 3 3
		t9 = 2 2 2 2 0 0 0 0 	*/
getelem(t5, ct, oph);
getelem(t2, ct, v2);   
if (iphase == 2)
	decr(v2);	/* only support hypercomplex FAD */

if (ix == 1)
	d2_init = d2;
t1_counter = (int) ((d2 - d2_init)*sw1 + 0.5);
if (t1_counter % 2)
	{
		add(v2, two, v2);
		add(oph, two, oph);
	}

status(A);
rlpower(pwxlvl, DODEV);
rcvroff();
if (sstrim != 0.0)
	{
		rlpower((tpwr - 4.0), TODEV);
		rgpulse(sstrim,zero,2.0e-5,rof1);
		rgpulse(sstrim,one,0.0,0.0);
	}
rlpower(tpwr, TODEV);
hsdelay(d1); 
status(B);
if (null != 0.0)
	{
		txphase(t1);
		rgpulse(pw, t1, rof1, 0.0);
		delay(bird - 2*pwx);
		decpulse(pwx, t8);
		simpulse(2*pw, 2*pwx, t1, t1, 2.0e-6, 2.0e-6);
		decpulse(pwx, t8);
		txphase(t9);
		delay(bird - 2*pwx);
		rgpulse(pw, t9, 0.0, rof1);
		hsdelay(null);
	}
txphase(t1);
rgpulse(pw, t1, rof1, 2.0e-6);
delay(bird);
simpulse(2*pw, pwx, t1, v2, 2.0e-6, 2.0e-6);
txphase(t3);
delay(bird);
simpulse(mult*pw, 2*pwx, t3, zero, 2.0e-6, 2.0e-6);
txphase(t1);
delay(bird);
delay(d2/2.0);
rgpulse(2*pw, t1, rof1, rof1);
decphase(t4);
delay(d2/2.0);
decpulse(pwx, t4);
rcvron(); 
rlpower(dpwr, DODEV);
delay(bird - 4.2e-6);
status(C);
}

