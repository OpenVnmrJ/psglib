
/* wlexch - solid-state echo sequence for exchange
  ************************************************
  method published by Bluemich & Spiess
	phase table required is <wlexch>
*/

#include <standard.h>

pulsesequence()
{
   double sw1,tau1,t1evol,tau2,magic,mix,phase;

/* initialize parameters */
   magic=getval("magic");
   mix=getval("mix");
   sw1=getval("sw1");
   phase= getval("phase");
   tau1 = getval("tau1");
   tau2 = getval("tau2");
   t1evol=d2;
   if (t1evol > 0.5/sw1)  t1evol=t1evol-(4.0*pw/3.1416);

	loadtable("wlexch");

   status(A);
      delay(d1);

   status(B);
      rcvroff();
      rgpulse(pw,t1, rof1, 0.0);
      delay(t1evol);
      if (phase == 1.0)	rgpulse(magic,t2,0.0,0.0);	   /*cosine*/
	else rgpulse(magic,t6,0.0,0.0);			   /*sine*/
      rcvron();
      delay(mix-rof1);
      if (phase == 1.0) rgpulse(magic,t3,rof1,0.0);	   /*cosine*/
	else  rgpulse(magic,t7,rof1,0.0);		   /*sine*/
      rgpulse(pw, t4, tau1, rof2);
      delay(tau2 - rof2 + pw);
/* if you have an older system then the status command is OK
   status(C);
but if the system has a type 3 AP Interface leave it commented out
as it will use 6usec otherwise */
      setreceiver(t5);
      delay(1.0 / (beta * fb));	/* wait out group delay */
      acquire(np, 1.0 / sw);	/* acquire data */
}
