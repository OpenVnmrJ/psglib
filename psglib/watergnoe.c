#ifndef LINT
static char     SCCSid[] = "watergnoe.c RC Burroughs Wellcome Co.";
#endif

/* watergnoe -  2D cross relaxation experiment.  It can be performed in
            either phase-sensitive or absolute value mode.  Either
            TPPI or the hypercomplex method can be used to achieve
            F1 quadrature in a phase-sensitive presentation.  No
            attempt is made to suppress J-cross peaks in this pulse
            sequence.  A watergate read pulse is employed.


  Parameters:

       mix = mixing time.
     phase =   0: gives non-phase-sensitive experiment (P-type peaks);
                  nt(min) = multiple of 16
                  nt(max) = multiple of 64

             1,2: gives HYPERCOMPLEX phase-sensitive experiment;
               3: gives TPPI phase-sensitive experiment;
                  nt(min) = multiple of  8
                  nt(max) = multiple of 32

	gtm	gradient duration during mix time
	glvlm	level for gtm zgradpulse   
	gt1	gradient duration during watergate read pulse (0.002 sec)
	glvl1	level for gt1 zgradpulse  (11000)	
	psel	selective 90 duration  (6090 usec)
	psellvl	power level for psel pulse (-2)
	phaseinc	0.5 degree phaseshift for softpulses  (1)

	if 13C channel and 15N channel are configured properly with
	powers modulation freqs, and offsets you can decouple labels
	with dm='nyny' and or dm2='nyny'.

	Macro defaults to dm & dm2 = 'nnnn'

	THIS SEQUENCE DOES MILD CHECKING OF dpwr and/or dpwr2 VALUES!
	IT IS YOUR RESPONSIBILITY TO NOT MELT THE PROBE.

 watergate noesy.c  31 Oct 94 R. Crouch
			crouch@bwco.com		 */


#include <standard.h>

pulsesequence()
{
   double          arraydim,
                   ss,sstrim,
                   psel,psellvl,gt1,glvl1,grise,
                   ni,gtm,glvlm,phaseinc,
                   mix;
   int             phase;


/* LOAD VARIABLES */
   ni = getval("ni");
   arraydim = getval("arraydim");
   mix = getval("mix");
   phase = (int) (getval("phase") + 0.5);
   ss = getval("ss");
   sstrim = getval("sstrim");
   psel = getval("psel");
   psellvl = getval("psellvl");
   gtm = getval("gtm");
   glvlm = getval("glvlm");
   gt1 = getval("gt1");
   glvl1 = getval("glvl1");
   grise = getval("grise");
  phaseinc = getval("phaseinc");
 
  if (phaseinc < 0.0) 
    phaseinc = 720+phaseinc;
  stepsize(0.5, TODEV);
  initval(phaseinc, v7);


/* For TPPI,  initialize v14 with the t1 increment number */
   if (phase == 3)
      initval((double) ((int) ((ix - 1) / (arraydim / ni) + 1e-6)), v14);


/* CHECK CONDITIONS */
   if ((rof1 < 9.9e-6) && (ix == 1))
      fprintf(stdout,"Warning:  ROF1 is less than 10 us\n");

/* Check for correct DM settings */
   if ((dm[A] == 'y') || (dm2[A] == 'y'))
   {
      printf("DM/DM2 must be set to either 'nyny', 'nnny', or 'nnnn'.\n");
      abort(1);
   }

/* Check dpwr & dpwr2 */
   if ((dpwr > 52) || (dpwr2 > 50))
   {
	printf("DPWR and/or DPWR2 out of range. Don't melt the probe.\n");
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


/* STEADY-STATE PHASECYCLING
/* This section determines if the phase calculations trigger off of (SS - SSCTR)
   or off of CT */

   ifzero(ssctr);
      mod2(ct, v2);
      hlv(ct, v3);
   elsenz(ssctr);
      sub(ssval, ssctr, v12);	/* v12 = 0,...,ss-1 */
      mod2(v12, v2);
      hlv(v12, v3);
   endif(ssctr);


/* CALCULATE PHASECYCLE */
   dbl(v2, v2);
   hlv(v3, v10);
   hlv(v10, v10);
   if (phase == 0)
   {
      assign(v10, v9);
      hlv(v10, v10);
      mod2(v9, v9);
   }
   else
   {
      assign(zero, v9);
   }
   hlv(v10, v1);
   mod2(v1, v1);
   dbl(v1, v1);
   add(v9, v2, v2);
   mod2(v10, v10);
   add(v1, v2, oph);
   add(v3, oph, oph);
   add(v10, oph, oph);
   add(v10, v1, v1);
   add(v10, v2, v2);
   add(v10, v3, v3);
   add(v3, two, v4);  /* for selective 90s in watergate read */

   if (phase == 2)
      incr(v2);
   if (phase == 3)
      add(v2, v14, v2);		/* TPPI phase increment */

/* The first 90 degree pulse is cycled first to suppress axial peaks.
   This requires a 2-step phasecycle consisting of (0 2).  The
   third 90 degree pulse is cycled next using a 4-step phasecycle
   designed to select both longitudinal magnetization, J-ordered states,
   and zero-quantum coherence (ZQC) during the mixing period.  If the
   experiment is to collect data requiring an absolute value display,
   the first pulse and the receiver are next incremented by 1 to achieve
   w1 quadrature (P-type peaks).  If the data are to be presented in
   a phase-sensitive manner, this step is not done.  Next, the second
   90 degree pulse is cycled to suppress axial peaks.  Finally, all
   pulse and receiver phases are incremented by 90 degrees to achieve
   quadrature image suppression due to receiver channel imbalance. */



/* BEGIN THE ACTUAL PULSE SEQUENCE */
   status(A);
   if (sstrim != 0.0)
   {
	rlpower((tpwr - 3.0), TODEV);
	rgpulse(sstrim, zero, 6.0e-5,2.0e-6);
	rgpulse(sstrim, one, 2.0e-6, rof1);
        rlpower(tpwr, TODEV);
	zgradpulse(glvl1, gt1);
   }
   delay(d1);
   rcvroff(); delay(rof1);
     status(B);
      rgpulse(pw, v2, rof1, 1.0e-6);
      if (d2 > 0.0)
      delay(d2 - 1.0e-6 - rof1 - (4.0*pw/3.14159));
      rgpulse(pw, v1, rof1, 1.0e-6);
   status(C);
      delay(mix - gtm - grise);
      zgradpulse(glvlm, gtm);
      delay(grise);
  /* watergate read pulse */
  pulse(pw,v3); delay(rof1);
 delay(grise);
 zgradpulse(glvl1, gt1);
 rlpower(psellvl, TODEV);
 xmtrphase(v7);
 delay(grise);
 txphase(v4);
 pulse(psel,v4);
 xmtrphase(zero);
 txphase(v3); rlpower(tpwr, TODEV);
 delay(2.0e-6);
 pulse(2*pw,v3);
 txphase(v4);
 rlpower(psellvl, TODEV); 
 xmtrphase(v7);
 delay(2.0e-6);
 pulse(psel,v4);
 delay(grise);
 zgradpulse(glvl1, gt1);
 rlpower(tpwr, TODEV);
 xmtrphase(zero);
 delay(grise);
 status(D);
 rcvron(); delay(rof1);
}

