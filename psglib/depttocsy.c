/* depttocsy.c  RC version 1.0 depthmqc followed by tocsy
	linear amps required */

#include <standard.h>

mleva()
{
	rgpulse(p1,v11,0.0,0.0);  delay(getval("window"));
	rgpulse(2*p1,v12,0.0,0.0); delay(getval("window"));
	rgpulse(p1,v11,0.0,0.0);
}

mlevb()
{
	rgpulse(p1,v13,0.0,0.0);  delay(getval("window"));
	rgpulse(2*p1,v10,0.0,0.0); delay(getval("window"));
	rgpulse(p1,v13,0.0,0.0);
}

	
static double	d2_init = 0.0;

pulsesequence()
{
	double	pwx, pwxlvl, j, multh, multc, bird, trim,
		p1lvl, mix, window, cycles, null, sw1, sstrim;
	int		iphase, t1_counter;

pwx = getval("pwx");
pwxlvl = getval("pwxlvl");
j = getval("j");
multh = getval("multh");
multc = getval("multc");
null = getval("null");
sw1 = getval("sw1");
p1lvl = getval("p1lvl");
sstrim = getval("sstrim");
mix = getval("mix");
trim = getval("trim");
window = getval("window");
iphase = (int)(getval("phase") + 0.5);

if (j != 0.0)
	bird = 1.0/(2.0*j);
else
	bird = 3.1e-3;

loadtable("depttocsy");
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
/* calculate mlev phases */
mod2(ct,v1);	/*v1 = 01 ..*/
dbl(v1,v1);	/*v1 = 02 ..*/
hlv(v1,v1);
hlv(v1,v1);	/*v1 = [0]4 [2]4 */
hlv(v1,v3);	/*v3 = [0]8 [2]8 */
add(v1,v3,v11);
add(v11,one,v12);
add(v12,one,v13);
add(v13,one,v10);

status(A);
rlpower(pwxlvl, DODEV);
rcvroff();
if (sstrim != 0.0)
	{
		rlpower((tpwr - 4.0), TODEV);
		rgpulse(sstrim,zero,5.0e-5,rof1);
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
		simpulse(2*pw, 2*pwx, t1, t1, 0.0, 0.0);
		decpulse(pwx, t8);
		txphase(t9);
		delay(bird - 2*pwx);
		rgpulse(pw, t9, 0.0, rof1);
		hsdelay(null);
	}
txphase(t1); delay(rof1);
rgpulse(pw, t1, rof1, 2.0e-6);
delay(bird);
simpulse(2*pw, pwx, t1, v2, 2.0e-6, 2.0e-6);
txphase(t3);
delay(bird);
simpulse(multh*pw, 2*pwx, t3, zero, 2.0e-6, 2.0e-6);
txphase(t1);
delay(bird);
delay(d2/2.0);
rgpulse(2*pw, t1, rof1, rof1);
decphase(t4);
delay(d2/2.0);
decpulse(pwx, t4);
rlpower(p1lvl, TODEV);
txphase(v12); decphase(zero);
rlpower(dpwr, DODEV);
delay(bird - 8.4e-6);
status(C);
cycles = (mix-trim)/(64.66*p1+32*window);
cycles = 2.0*(double)(int)(cycles/2.0);
initval(cycles,v5);
rgpulse(trim,v12,5.0e-6,0.0);
if (cycles > 1.0)
	{
	 starthardloop(v5);
		mleva(); mlevb(); mlevb(); mleva();
		mlevb(); mlevb(); mleva(); mleva();
		mlevb(); mleva(); mleva(); mlevb();
		mleva(); mleva(); mlevb(); mlevb();
		rgpulse(0.66*p1,v12,0.0,0.0);
	 endhardloop();
	}
if (multc != 0.0)
	{		/* suppression or inversion of direct responses */
		status(D);
		rlpower(pwxlvl, DODEV);
		rlpower(tpwr, TODEV);
		delay(bird - 8.4e-6 - rof1);
		simpulse(2*pw,multc*pwx,t8,t1,rof1,rof1);
		delay(bird - 4.2e-6);
	}
rlpower(dpwr, DODEV); decphase(zero); rcvron(); delay(rof1);
status(E);
}
