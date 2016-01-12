/* slpdqcosy - double quantum filtered cosy experiment
              with SLP presaturation option.

    [Added features : 1. FAD-hypercomplex
                      2. xmtr presaturation       ]


   ref: u. piantini, o.w. sorenson, and r.r. ernst,
        j. am. chem. soc. 104:6800-6801 (1982)
        m. rance et al., bbrc 117:479-485 (1983)


  Parameters:

         pw = 90 excitation pulse (at power level tpwr)
      phase =   0: P-type non-phase-sensitive experiment
              1,2: hypercomplex phase-sensitive experiment
                3: TPPI phase-sensitive experiment
     satflg = 'y' uses obs xmtr for presat at satfrq and satpwr for 
              satdly seconds
     satfrq = saturation frequency for xmtr presaturation
              (used only if different from tof)
     satpwr = saturation power for xmtr presaturation
     satdly = saturation period follows D1
     slpflg = flag that turns on slp presaturation
    satshape= shape used for slp presaturation
               slpsatd (tof > satfrq)
               slpsatu (tof < satfrq)
         ss = number of steady state pulses; if ss < 0 then -ss
              steady-state transients are performed before EVERY increment;
              if ss > 0, then ss steady-state transients are performed
              before only the first increment.
      sspul = 'y': selects for trim(x)-trim(y) at start of pulse sequence
              'n': normal D1 delay

         nt = min:  multiple of 16 (phase=0)
                    multiple of 8  (phase=1,2  phase=3)
              max:  multiple of 64 (phase=0)
                    multiple of 32 (phase=1,2  phase=3)


  NOTE:  If phase = 3, remember that sw1 must be set to twice the
         desired value.  The 4-january revision included the following
         sequence at the beginning of the pulse sequence:  homospoil -
         90 degree pulse - homospoil.  This should eliminate both the
         DQ-like artifacts in the 2D spectrum and the oscillatory nature
         of the steady-state.  This inclusion is selected if sspul='y'.


  s. farmer     21 july       1987
  revised       21 december   1987
  revised        4 january    1988
  revised       21 january    1988 
  
  revised       07 May 1990 (vvk)
  revised       23 Apr 1992 (vvk) */


#include <standard.h>

pulsesequence()
{
   double          ss,
                   phase,
                   satdly,
                   satfrq,
                   cycle,
                   shape_pw,
                   h2off,
                   satpwr;
   int             iphase;
   char            sspul[MAXSTR],
                   slpflg[MAXSTR],
                   satshape[MAXSTR],
                   satflg[MAXSTR];


/* LOAD VARIABLES AND CHECK CONDITIONS */
   satdly = getval("satdly");
   satpwr = getval("satpwr");
   h2off = getval("h2off");
   satfrq = getval("satfrq");
   phase = getval("phase");
   iphase = (int) (phase + 0.5);
   getstr("sspul", sspul);
   getstr("satflg",satflg);
   getstr("slpflg",slpflg);
   getstr("satshape",satshape);

   if ((slpflg[0] == 'y') || (slpflg[1] == 'y'))
   {
    if (h2off != 0.0)
    {
     shape_pw = 10.0/h2off;
     if (shape_pw < 0.0)
      shape_pw = -shape_pw;

     cycle = (satdly/(shape_pw + 15.4e-6));
     cycle = 2.0*(double)(int)(cycle/2.0);
     initval(cycle,v10);
    }
    if (rfwg[0] != 'y')
    {
     printf("NOTE: slpflg=y requires PPM in Channel 1. psg aborted.\n");
     abort(1);
    }
   }

   if ((slpflg[1] == 'y') && (h2off != 0.0))
   {
    printf("NOTE: off-resonance presat cannot be done during d2.\n");
    abort(1);
   }

   if ((satflg[1] == 'y') && (satfrq != tof))
   {
    printf("NOTE: off-resonance presat cannot be done during d2.\n");
    abort(1);
   }

   ss = getval("ss");

   if (iphase == 3)
      initval((double) (ix - 1), v14);
   else
   {
    if ((iphase == 1) || (iphase == 2))
      initval(2.0*(double)(((int)(d2*getval("sw1")+0.5)%2)),v14);
    else
      assign(zero, v14);
   }

   if ((rof1 < 9.9e-6) && (ix== 1))
      fprintf(stdout,"Warning:  ROF1 is less than 10 us\n");



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

   if ((iphase == 1) || (iphase == 2))
     add(oph,v14,oph);


/* BEGIN ACTUAL PULSE SEQUENCE CODE */
   status(A);
   rlpower(tpwr,TODEV);
   delay(5.0e-6);
   if (sspul[0] == 'y')
   {
      rgpulse(200*pw, zero, rof1, 0.0);
      rgpulse(200*pw, one, 0.0, rof1);
   }

   hsdelay(d1);

   if (satflg[0] == 'y')
   {
      rlpower(satpwr,TODEV);
      if (slpflg[0] == 'y')
      {
       if (h2off != 0.0)
       {
        rcvroff();
        delay(rof1);
        starthardloop(v10);
         shaped_pulse(satshape,shape_pw,zero,0.0,1e-6);
        endhardloop();
        delay(rof2);
        rcvron();
       }
       else
       {
        if (satfrq != tof)
         offset(satfrq, TODEV);
        shaped_pulse(satshape,satdly,zero,rof1,rof2);
        if (satfrq!= tof)
         offset(tof, TODEV);
       }
      }
      else
      {
       if (satfrq != tof)
        offset(satfrq, TODEV);
       rgpulse(satdly,zero,rof1,rof2);
       if (satfrq!= tof)
        offset(tof, TODEV);
      }
      rlpower(tpwr,TODEV);
      delay(40e-6);
   }

   status(B);
      rgpulse(pw, v1, rof1, 1.0e-6);

      if (satflg[1] == 'y')
      {
       rlpower(satpwr,TODEV);
       if (slpflg[1] == 'y')
       {
        if (d2 > 0.0)
         shaped_pulse(satshape,d2 - rof1 - 24.8e-6 - (4*pw)/3.1414,zero,0.0,0.0);
       }
       else
       {
        if (d2 > 0.0)
         rgpulse(d2 - rof1 - 9.4e-6 - (4*pw)/3.1414,zero,0.0,0.0);
       }
       rlpower(tpwr,TODEV);
      }
      else
      if (d2 > 0.0)
       delay(d2 - rof1 - 1.0e-6 - (4*pw)/3.1414);

      rgpulse(pw, v2, rof1, 0.0);
      rgpulse(pw, v3, 1.0e-6, rof2);
   status(C);
}
