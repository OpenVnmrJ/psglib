#ifndef LINT
static char SCCSid[] = "@(#)dqcosy.c 9.1 9/28/93 Copyright (c) 1991,1993 Varian Assoc.,Inc. All Rights Reserved";
#endif

/* wetdqcosy.c - made from dqcosy.c (V9.1 9/28/93)
   dqcosy - double quantum filtered cosy experiment


   ref: u. piantini, o.w. sorenson, and r.r. ernst,
        j. am. chem. soc. 104:6800-6801 (1982)
        m. rance et al., bbrc 117:479-485 (1983)


  Parameters:

         pw = 90 excitation pulse (at power level tpwr)
      phase =   0: P-type non-phase-sensitive experiment
              1,2: hypercomplex phase-sensitive experiment
                3: TPPI phase-sensitive experiment
      sspul = 'y': selects for HS-90-HS sequence at start of pulse sequence
              'n': normal DQCOSY experiment

         nt = min:  multiple of 16 (phase=0)
                    multiple of 8  (phase=1,2  phase=3)
              max:  multiple of 64 (phase=0)
                    multiple of 32 (phase=1,2  phase=3)
	pwwet
	wetpwr
	wetshape
	gswet
	c13wet
	dmfwet
	dpwrwet
	dofwet


  NOTE:  If phase = 3, remember that sw1 must be set to twice the
         desired value.  The 4-january revision included the following
         sequence at the beginning of the pulse sequence:  homospoil -
         90 degree pulse - homospoil.  This should eliminate both the
         DQ-like artifacts in the 2D spectrum and the oscillatory nature
         of the steady-state.  This inclusion is selected if sspul='y'.


  P.A.Keifer 950407 - made standard format wetdqcosy.c
  P.A.Keifer 950920 - updated wet
  P.A.Keifer 960116 - added tpwrf control to wet4


	 */


#include <standard.h>

pulsesequence()
{
   double          phase,
                   corr;
   int             iphase;
   char            sspul[MAXSTR];


/* LOAD VARIABLES AND CHECK CONDITIONS */
   phase = getval("phase");
   iphase = (int) (phase + 0.5);
   getstr("sspul", sspul);

   if (phase > 2.5)
   {
      initval((double) (ix - 1), v14);
   }
   else
   {
      assign(zero, v14);
   }

   if ((rof1 < 9.9e-6) && (ix== 1))
      fprintf(stdout,"Warning:  ROF1 is less than 10 us\n");

/* STEADY-STATE PHASECYCLING
/* This section determines if the phase calculations trigger off of (SS - SSCTR)
   or off of CT */

   ifzero(ssctr);
      hlv(ct, v4);
      mod4(ct, v3);
   elsenz(ssctr);
      sub(ssval, ssctr, v12);	/* v12 = 0,...,ss-1 */
      hlv(v12, v4);
      mod4(v12, v3);
   endif(ssctr);


/* CALCULATE PHASECYCLE */
/* The phasecycle first performs a 4-step cycle on the third pulse in order
   to select for DQC.  Second, the 2-step QIS cycle is added in.  Third, a
   2-step cycle for axial peak suppression is performed on the second pulse.
   Fourth, a 2-step cycle for axial peak suppression is performed on the
   first pulse.  If P-type peaks only are being selected, the 2-step cycle
   for P-type peak selection is performed on the first pulse immediately
   after the 4-step cycle on the third pulse. */

   hlv(v4, v4);
   if (iphase == 0)
   {
      assign(v4, v6);
      hlv(v4, v4);
      mod2(v6, v6);		/* v6 = P-type peak selection in w1 */
   }
   hlv(v4, v2);
   mod4(v4, v4);		/* v4 = quadrature image suppression */
   hlv(v2, v1);
   mod2(v1, v1);
   dbl(v1, v1);
   mod2(v2, v2);
   dbl(v2, v2);
   dbl(v3, v5);
   add(v3, v5, v5);
   add(v1, v5, v5);
   add(v2, v5, v5);
   add(v4, v5, v5);
   add(v4, v1, v1);
   add(v4, v2, v2);
   add(v4, v3, v3);
   if (iphase == 0)
   {
      add(v6, v1, v1);
      add(v6, v5, v5);
   }
   if (iphase == 2)
      incr(v1);
   add(v14, v1, v1);		/* adds TPPI increment to the phase of the
				 * first pulse */
   assign(v5, oph);
/* Add FAD for phase=1 or phase=2 */
   if ((iphase == 1) || (iphase == 2))
   {
      initval(2.0*(double)((int)(d2*getval("sw1")+0.5)%2),v13);
      add(v1,v13,v1); add(oph,v13,oph);
   }


/* BEGIN ACTUAL PULSE SEQUENCE CODE */
   status(A);
   if (sspul[0] == 'y')
   {
      hsdelay(hst + 0.001);
      rgpulse(pw, v1, 1.0e-6, 1.0e-6);
      hsdelay(hst + 0.001);
   }
      hsdelay(d1);
    if (getflag("wet")) wet4(zero,one);
   status(B);
      rgpulse(pw, v1, rof1, 1.0e-6);
      corr = 1.0e-6 + rof1 + 4.0*pw/3.1416;
      if (d2  > corr)
        delay(d2-corr);
      rgpulse(pw, v2, rof1, 0.0);
      rgpulse(pw, v3, 1.0e-6, rof2);
   status(C);
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


