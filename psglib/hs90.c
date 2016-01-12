
/* flipflop - setup sequences for multipulse

	This sequence will run both flip flip and flip flop

	( pw(phase1)-dtau-acq-pw(phase2)-dtau-)np/4
          l------------------l
                   tau


   Contact: e.h.williams, palo alto mar 1989
   Version:  revision of flipflop.c

   tau is cycle time
   phase1 is phase of 1st pulse
   phase2 is phase of 2nd
   phaser is the phase of the Tx compared with the Rx.

   IPA is provided on tpwrf

 */


#include <standard.h>

pulsesequence()
{
   double          tau,
                   dtau,
                   phase1,
                   phase2,
		   phaser,
                   points;
   char            trig[MAXSTR];

/* IPA definition */
   ipwrf(tpwrf,TODEV,"tpwrf");

   stepsize(0.5, TODEV);
   tau = getval("tau");
   phaser = getval("phaser");
   phase1 = getval("phase1");
   phase2 = getval("phase2");
   getstr("trig", trig);

   points = 2.0;		/* 2 points per acq is normal */
   dtau = tau - pw - rof1 - rof2 - (1.0e-7) * points;

   initval(np / 2.0, v9);
   initval(phase1, v1);
   initval(phase2, v2);
   if (phaser < 0.0) {phaser=phaser+360.0;}     /*make shift positive */
   phaser= (double) ( (int) (phaser*2.0+0.9));	/* calculate phase shift*/
   initval(phaser,v8);
   stepsize(0.5,TODEV);

   if (d1 < 1.0)
   {
      fprintf(stdout, "d1 delay is too short!\n");
      abort(1);
   }

   status(A);
   xmtrphase(zero);
   xmtrphase(v8);	/*set small angle phase */
   delay(d1);


   if (trig[0] == 'y')
   {
      xgate(1.0);
      fprintf(stdout, "Synchronized acquisition.\n");
   }

   starthardloop(v9);

      rgpulse(pw, v1, rof1, rof2);
      delay(dtau);
      delay(0.2e-6);

      rgpulse(pw, v2, rof1, rof2);
      delay(dtau);
      delay(0.2e-6);

      rgpulse(pw, v2, rof1, rof2);
      delay(dtau);
      delay(0.2e-6);

      rgpulse(pw, v1, rof1, rof2);
      delay(dtau);
      delay(0.2e-6);

      rgpulse(pw, v1, rof1, rof2);
      delay(dtau);
      delay(0.2e-6);

      rgpulse(pw, v2, rof1, rof2);
      delay(dtau);
      delay(0.2e-6);

      rgpulse(pw, v2, rof1, rof2);
      delay(dtau);
      delay(0.2e-6);

      rgpulse(pw, v1, rof1, rof2);
      delay(dtau);
      acquire(points, 0.2e-6);

   endhardloop();
}
