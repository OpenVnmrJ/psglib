#ifndef LINT
static char SCCSid[] = "@(#)noesy.c 9.2 10/11/94 Copyright (c) 1991,1993 Varian Assoc.,Inc. All Rights Reserved";
#endif

/* wetnoesy.c - made from noesy.c (v9.2 10/11/94)
   noesy -  2D cross relaxation experiment.  It can be performed in
            either phase-sensitive or absolute value mode.  Either
            TPPI or the hypercomplex method can be used to achieve
            F1 quadrature in a phase-sensitive presentation.  No
            attempt is made to suppress J-cross peaks in this pulse
            sequence.


  Parameters:

       mix = mixing time.
     phase =   0: gives non-phase-sensitive experiment (P-type peaks);
                  nt(min) = multiple of 16
                  nt(max) = multiple of 64

             1,2: gives HYPERCOMPLEX phase-sensitive experiment;
               3: gives TPPI phase-sensitive experiment;
                  nt(min) = multiple of  8
                  nt(max) = multiple of 32

    presat = length of decoupler presaturation period; if PRESAT > 0,
             D1 is reduced to (D1 - PRESAT) and the presaturaion period
             is added in after that period; does not depend on DM but
             does depend on DMM.
		DISABLED HERE !!!
     sspul = 'y': selects for HS-90-HS sequence at start of pulse sequence
             'n': normal NOESY experiment


	P.A.Keifer 950406 - made wetnoesy.c
	P.A.Keifer 950920 - updated wet
        P.A.Keifer 960116 - added tpwrf control to wet4 and added 45 degree fix

 */


#include <standard.h>

pulsesequence()
{
   double          arraydim,
                   corr,
                   mix,pwwet,gtw,gswet,dz;
   int             phase;
   char            sspul[MAXSTR];


/* LOAD VARIABLES */
  pwwet=getval("pwwet");        /* User enters power for 90 deg. */
  gtw=getval("gtw");            /* Z-Gradient duration           */
  gswet=getval("gswet");        /* Post-gradient stability delay */
  dz=getval("dz");
   ni = getval("ni");
   arraydim = getval("arraydim");
   mix = getval("mix");
   phase = (int) (getval("phase") + 0.5);
   getstr("sspul", sspul);

/* For TPPI,  initialize v14 with the t1 increment number */
   if (phase == 3)
      initval((double) ((int) ((ix - 1) / (arraydim / ni) + 1e-6)), v14);


/* CHECK CONDITIONS */
   if ((rof1 < 9.9e-6) && (ix == 1))
      fprintf(stdout,"Warning:  ROF1 is less than 10 us\n");
   if ((mix - rof1 - (4.0*pwwet) - (4.0*20.0e-6) - (4.0*rof2) - (4.0*gtw) - (4.0*gswet)) < 1.0e-6)
      {
      printf("Warning: mix time is too short for WET to be executed.\n");
      abort(1);
      }

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
   if (phase == 2)
      incr(v2);
   if (phase == 3)
      add(v2, v14, v2);		/* TPPI phase increment */

/* FAD added for phase=1 or phase=2 */
   if ((phase == 1) || (phase == 2))
   {
      initval(2.0*(double)((int)(d2*getval("sw1")+0.5)%2),v11);
      add(v2,v11,v2); add(oph,v11,oph);
   }

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
     initval(7.0,v7);
     stepsize(45.0,TODEV);
     xmtrphase(v7);
     delay(0.2e-6);
   if (sspul[0] == 'y')
   {
      hsdelay(hst + 0.001);
      rgpulse(pw, v2, 1.0e-6, 1.0e-6);
      hsdelay(hst + 0.001);
   }
      hsdelay(d1);

   status(B);
      rgpulse(pw, v2, rof1, 1.0e-6);
      xmtrphase(zero);
      delay(0.2e-6);
      corr = 0.2e-6 + 1.0e-6 + rof1 + 4.0*pw/3.1416;
      if (d2  > corr)
        delay(d2-corr);
      rgpulse(pw, v1, rof1, 1.0e-6);
   status(C);
    if (getflag("wet"))
	{
          hsdelay(mix - rof1 - (4.0*pwwet) - (4.0*20.0e-6) - (4.0*rof2) - (4.0*gtw) - (4.0*gswet) - dz);
	 wet4(zero,one);
	}
	else
          hsdelay(mix - rof1);

   status(D);
      rgpulse(pw, v3, rof1, rof2);
}
/* wet4 - Water Elimination */
wet4(phaseA,phaseB)
  codeint phaseA,phaseB;

{
  double finepwr,gzlvlw,gtw,gswet,dmfwet,dpwrwet,dofwet,wetpwr,pwwet,dz;
  int c13wet;
  char   wetshape[MAXSTR];
  c13wet=getflag("c13wet");             /* Water suppression flag        */
  getstr("wetshape",wetshape);    /* Selective pulse shape (base)  */
  wetpwr=getval("wetpwr");        /* User enters power for 90 deg. */
  pwwet=getval("pwwet");        /* User enters power for 90 deg. */
  dmfwet=getval("dmfwet");
  dpwrwet=getval("dpwrwet");
  dofwet=getval("dofwet");
  dz=getval("dz");
  finepwr=wetpwr-(int)wetpwr;     /* Adjust power to 152 deg. pulse*/
  wetpwr=(double)((int)wetpwr);
  if (finepwr==0.0) {wetpwr=wetpwr+5; finepwr=4095.0; }
  else {wetpwr=wetpwr+6; finepwr=4095.0*(1-((1.0-finepwr)*0.12)); }
  rcvroff();
  if (c13wet)
    {
    setstatus(DECch,FALSE,'w',FALSE,dmfwet);
    decoffset(dofwet);
    decpower(dpwrwet);
    }
  obspower(wetpwr);         /* Set to low power level        */
  gzlvlw=getval("gzlvlw");      /* Z-Gradient level              */
  gtw=getval("gtw");            /* Z-Gradient duration           */
  gswet=getval("gswet");        /* Post-gradient stability delay */
  chess(finepwr*0.5059,wetshape,pwwet,phaseA,20.0e-6,rof2,gzlvlw,gtw,gswet,c13wet);
  chess(finepwr*0.6298,wetshape,pwwet,phaseB,20.0e-6,rof2,gzlvlw/2.0,gtw,gswet,c13wet);
  chess(finepwr*0.4304,wetshape,pwwet,phaseB,20.0e-6,rof2,gzlvlw/4.0,gtw,gswet,c13wet);
  chess(finepwr*1.00,wetshape,pwwet,phaseB,20.0e-6,rof2,gzlvlw/8.0,gtw,gswet,c13wet);
  if (c13wet)
    {
    setstatus(DECch,FALSE,'c',FALSE,dmf);
    decoffset(dof);
    decpower(dpwr);
    }
  obspower(tpwr);               /* Reset to normal power level   */
  obspwrf(tpwrf);
  rcvron();
  delay(dz);
}
 
/* chess - CHEmical Shift Selective Suppression */
chess(pulsepower,pulseshape,duration,phase,rx1,rx2,gzlvlw,gtw,gswet,c13wet)  double pulsepower,duration,rx1,rx2,gzlvlw,gtw,gswet;
  int c13wet;
  codeint phase;
  char* pulseshape;
{
  obspwrf(pulsepower);
  if (c13wet) decon();
  shaped_pulse(pulseshape,duration,phase,rx1,rx2);
  if (c13wet) decoff();
  zgradpulse(gzlvlw,gtw);

  delay(gswet);
}

int getflag(str)
char str[MAXSTR];
{
   char strval[MAXSTR];

   getstr(str,strval);
   if ((strval[0]=='y') || (strval[0]=='Y')) return(TRUE);
     else                                    return(FALSE);
}

 
