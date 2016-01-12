/* BIRDHMBC - Gradient Selected absolute value HMBC

        BIRDHMQC (gradient HMBC with excellent suppression of
	residual 1-bond correlations, superior to gHMBC, but 5 - 10%
	less sensitive than the latter), as published by Peter Bigler
	in Magn. Reson. Chem. 38, 963-969 (2000).

	Features included:
		F1 Axial Displacement
		Randomization of Magnetization prior to relaxation delay
			with G-90-G 
			[selected by sspul flag]
		J-filter to suppress one-bond correlations
				
	Paramters:
		sspul :		y - selects magnetization randomization option
		hsglvl:		Homospoil gradient level (DAC units)
		hsgt	:	Homospoil gradient time
		gzlvl1	:	encoding Gradient level
		gt1	:	encoding gradient time
		gzlvl3	:	decoding Gradient level
		gt3	:	decoding gradient time
		gstab	:	recovery delay
		j1xh	:	One-bond XH coupling constant
		jnxh	:	multiple bond XH coupling constant
		pwxlvl  :	X-nucleus pulse power
		pwx	:	X-nucleus 90 deg pulse width
		d1	:	relaxation delay
		d2	:	Evolution delay

Bruce Adams, 2001-02-26

*/


#include <standard.h>

static int ph3[4] = {1,1,3,3};
static int ph4[8] = {0,0,0,0,2,2,2,2};
static int ph6[8] = {0,2,0,2,2,0,2,0};

pulsesequence()
{
  double j1xh,
         jnxh,
	 pwxlvl,
	 pwx,
	 gzlvl1,
	 gt1,
	 gzlvl3,
	 gt3,
	 gstab,
	 hsglvl,
	 hsgt,
	 tau,
         taumb;
  char	 sspul[MAXSTR];

  j1xh = getval("j1xh");
  jnxh = getval("jnxh");
  pwxlvl = getval("pwxlvl");
  pwx = getval("pwx");
  getstr("sspul",sspul);
  gzlvl1 = getval("gzlvl1");
  gt1 = getval("gt1");
  gzlvl3 = getval("gzlvl3");
  gt3 = getval("gt3");
  gstab = getval("gstab");
  hsglvl = getval("hsglvl");
  hsgt = getval("hsgt");
  taumb = 1/(4.0*jnxh);
  tau=1/(2.0*j1xh);

  settable(t3,4,ph3);
  settable(t4,8,ph4);
  settable(t6,8,ph6);
 
  setreceiver(t6);

  /* initval(2.0*(double)((int)(d2*getval("sw1")+0.5)%2),v10);
  add(v3,v10,v3);
  add(oph,v10,oph); */

  mod2(ct,v1);

  status(A);
     decpower(pwxlvl);
     if (sspul[0] == 'y')
     {
	zgradpulse(hsglvl,hsgt);
	rgpulse(pw,zero,rof1,rof1);
	zgradpulse(hsglvl,hsgt);
     }

     delay(d1);
     /* rcvroff(); */

  status(B);
     rgpulse(pw,zero,rof1,rof1);
     delay(taumb - 2.0*rof1); 
     rgpulse(pw,zero,rof1,rof1);
     delay(tau - 2.0*rof1); 
     simpulse(2.0*pw, 2.0*pwx, t3, zero, rof1, rof1);
     delay(tau - 2.0*rof1); 
     rgpulse(pw,zero,rof1,rof1);
     delay(taumb - 2.0*rof1); 

     ifzero(v1);
       decrgpulse(2.0*pwx, zero, rof1, rof1);
     elsenz(v1);
       delay(2.0*(pwx+rof1));
     endif(v1);

     delay(tau - 2.0*rof1); 
     decrgpulse(pwx,t4,rof1,rof1);
     delay(d2/2.0);
     zgradpulse(gzlvl1,gt1);
     delay(gstab);
     rgpulse(2.0*pw, zero, rof1, rof1);
     zgradpulse(gzlvl1,gt1);
     delay(gstab);
     delay(d2/2);
     decrgpulse(pwx,zero,rof1,rof1);
     zgradpulse(gzlvl3,gt3);
     decpower(dpwr);
     /* rcvron(); */
     delay(gstab);
 
  status(C);
} 

