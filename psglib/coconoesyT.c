/* COCONOESYT.c  concurrent COSY & NOESY  in av mode.
	From J. Magn. Reson. 56, 343-349 (1984)
        Also Magnetic Moments Vol. II, No.3 p. 1 1986.
	See Magnetic Moments Vol. V, No. 2, p 6 (1992) for
	a description of the tablib-based random mix delay.

	Version 1.3 R.C. Crouch Burrough Wellcome Co.
	(919) 248-3840 or crouch@bwco.com              */

#include <standard.h>


pulsesequence()
{
  double           ss,
                   at,mixvar,
                   mix,sstrim;
 


/* LOAD VARIABLES */
   mix = getval("mix");
   mixvar = getval("mixvar");
   ss = getval("ss");
   fb = getval("fb");
   sw = getval("sw");
   nf = getval("nf");
   at = getval("at");
   sstrim = getval("sstrim");


/* CHECK CONDITIONS */
   if ((rof1 < 9.9e-6) && (ix == 1))
      fprintf(stdout,"Warning:  ROF1 is less than 10 us\n");

   if (at > mix)
        {
         printf("at too large  - acquisition aborted./n");
         abort(1);
        }

/* DETERMINE STEADY-STATE MODE */
   if (ss < 0)
   {
      ss = (-1) * ss;
   }
   else
   {
      if ((ss > 0) && (ix == 1))
      {
	 ss = ss;
      }
      else
      {
	 ss = 0;
      }
   }
   initval(ss, ssctr);
   initval(ss, ssval);

    if (mixvar != 0.0)
	 loadtable("coconoesyT");
    else
	loadtable("coconoe");

	/*
		t1 = 0 0 0 0 1 1 1 1 2 2 2 2 3 3 3 3
		t2 = 0 2 1 3 1 3 2 0 2 0 3 1 3 1 0 2
		t4 = 0 0 2 2 1 1 3 3 2 2 0 0 3 3 1 1
		t10 is necessary for random mixtime control
		e.g.
		t10 = 100 will give a fixed delay = mix.
		t10 = { 95 100 101 97 104 102 96 99 100
			103 98 96 97 104 96 98 101 99
			102 96 99 103 104 96 102 100 96
			98 100 99 97 103 95 104 98 101
			96 103 99 100 97}16
		would result in a random variation.
	*/

/* SETUP for the tablib-generated random mix time */
 
   initval(nt*d2*getval("sw1"), v8);
   add(v8,ct,v13);
   getelem(t10,v13,v11);

/* BEGIN THE ACTUAL PULSE SEQUENCE */
   status(A);
      power(zero,DODEV);
      if (sstrim != 0.0)
       {
        rlpower((tpwr - 4.0), TODEV);
        rgpulse(sstrim, zero, rof1, 0.0e-6);
        rgpulse(sstrim, one, 0.0e-6, rof1);
       }
      rlpower(tpwr, TODEV);
      setreceiver(t4);
      if (d1>hst) hsdelay(d1);
   status(B);
      rgpulse(pw, t1, rof1, 1.0e-6);
         delay(d2);
      rgpulse(pw, t2, 1.0e-6, 1.0e-6);
   status(C);
      if (nf == 2)
	{
          delay(alfa + 1.0/(2.0*fb));
	  acquire(np, 1.0/sw); 		     /* acquire cosy  */
	  delay(1.0e-2);	/* provide a little time for Mr. fifo */
          starthardloop(v11);
		delay(mix*0.005);
		delay(mix*0.005);
	  endhardloop();
	}
      else
	{
          starthardloop(v11);
		delay(mix*0.005);
		delay(mix*0.005);
	  endhardloop();
	}
      rgpulse(pw, t2, rof1, rof2);
    status(D);
    delay(alfa + 1.0/(2.0*fb));
    acquire(np, 1.0/sw);		/* acquire noesy */
}
