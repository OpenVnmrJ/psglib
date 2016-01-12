/* COCONOESY.c  concurrent COSY & NOESY  in av mode.
	From J. Magn. Reson. 56, 343-349 (1984)
        Also Magnetic Moments Vol. II, No.3 p. 1 1986.

	Version 1.2 R. Crouch Burrough Wellcome Co.
	(919) 248-3840 or crouch@bwco.com              */

#include <standard.h>

/* randomdelay-generate a random delay
  for different increments. Usage:
    randomdelay(time,variation), for
  example: randomdelay(mix,.1);
  where variation is a fraction, i.e. 10% */

#define constant 1073741824

randomdelay(time,variation)
 double time,variation;
{
  double value;
/* initialize on first increment */
 if (ix==1) srandom(getpid());
 value=(random()-constant)/(1.0*constant);
 hsdelay(time*(1.0+value*variation));
}


pulsesequence()
{
  double           ss,
                   at,
                   mix,mixvar,sstrim;
 


/* LOAD VARIABLES */
   mix = getval("mix");
   mixvar = getval("mixvar");   /* variation for random delay  */
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

	 loadtable("coconoesy");
	/*
		t1 = 0 0 0 0 1 1 1 1 2 2 2 2 3 3 3 3
		t2 = 0 2 1 3 1 3 2 0 2 0 3 1 3 1 0 2
		t4 = 0 0 2 2 1 1 3 3 2 2 0 0 3 3 1 1
	*/
 

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
          randomdelay((mix - at), mixvar);
	}
      else
       	 randomdelay(mix,mixvar);
      rgpulse(pw, t2, rof1, rof2);
    status(D);
    delay(alfa + 1.0/(2.0*fb));
    acquire(np, 1.0/sw);		/* acquire noesy */
}
