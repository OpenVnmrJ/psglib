
/* jsdqcosy -j-scaled partially-decoupled dqcosy experiment
               using obs xmtr for saturation
              TRANSMITTER MUST BY AT SOLVENT FREQUENCY

  Parameters:

      sspul = 'y': selects for Trim(x)-Trim(y)
               sequence at start of pulse sequence
             (highly recommended to eliminate "long T1" artifacts)
         pw = 90 excitation pulse (at power level tpwr)
      selpw = pulse width of selective inversion pulse
     selpwr = power level for selective inversion pulse
   selshape = shapelib file for selective inversion pulse 
               (must be in shapelib as .RF file. Should be created
                via "convolute" macro to make an "SLP" pulse, allowing
                off-resonance excitation without moving frequency)
   satshape = shapelib file for saturation (can use less power than with
              just attenuators)
    satmode = 'yn' for saturation during relaxation period
              'ny' for saturation during d2
              'yy' for both
     satdly = presaturation period -uses obs xmtr at tof with satpwr; 
      phase =   0: P-type non-phase-sensitive experiment
              1,2: hypercomplex phase-sensitive experiment
                3: TPPI phase-sensitive experiment
 

  revised from tndqcosy.c       Sept 1991  g.gray  */


#include <standard.h>

pulsesequence()
{
   double          ss,selpw,selpwr,
                   phase,
                   satdly,satpwr;
   int             iphase;
   char      satmode[MAXSTR], satshape[MAXSTR], selshape[MAXSTR], sspul[MAXSTR];


/* LOAD VARIABLES AND CHECK CONDITIONS */
   satdly = getval("satdly"); selpw=getval("selpw");
   phase = getval("phase"); selpwr=getval("selpwr");
   iphase = (int) (phase + 0.5);
   getstr("sspul", sspul); getstr("selshape",selshape);
   getstr("satmode", satmode); getstr("satshape",satshape);
   tpwr=getval("tpwr"); initval(tpwr,v10); initval(selpwr,v13);
   satpwr=getval("satpwr"); initval(satpwr,v11);
   ss = getval("ss");

       if (satpwr > 40)
        {
         printf("satpwr too large  - acquisition aborted./n");
         abort(1);
        }

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

  /* FOR HYPERCOMPLEX, USE REDFIED TRICK TO MOVE AXIALS TO EDGE */  
   initval(2.0*(double)(((int)(d2*getval("sw1")+0.5)%2)),v9);  /* moves axials */
   if ((iphase==2)||(iphase==1)) {add(v1,v9,v1); add(oph,v9,oph);}

/* BEGIN ACTUAL PULSE SEQUENCE CODE */
   status(A);
   if (sspul[0] == 'y')
   {
      rgpulse(200*pw, one, 10.0e-6, 0.0e-6);
      rgpulse(200*pw, zero, 0.0e-6, 1.0e-6);
   }
      hsdelay(d1);
    if (satmode[A] == 'y')
     {
      power(v11,TODEV); 
      shaped_pulse(satshape,satdly,zero,rof1,rof2);
      power(v10,TODEV);
     }
   status(B);
      rgpulse(pw, v1, rof1, 1.0e-6);
      if (satmode[B] == 'y')
       {
        power(v11,TODEV);
        if ((d2/2)>0)
           shaped_pulse(satshape,d2/2 -6e-6 -rof1 -(2*pw/3.1416),zero,5e-6,0.0);
        power(v13,TODEV);
        shaped_pulse(selshape,selpw,zero,rof1,rof1);
        power(v11,TODEV);
        if ((d2/4)>0)
           shaped_pulse(satshape,d2/4 -2*rof1 -2*pw -5e-6,zero,5e-6,rof1);
        power(v10,TODEV);
        rgpulse(pw,zero,rof1,0.0);
        rgpulse(2*pw,one,0.0,0.0);
        rgpulse(pw,zero,0.0,rof1);
        power(v11,TODEV);
        if ((d2/4)>0)
           shaped_pulse(satshape,d2/4 -2*rof1 -5e-6 -2*pw -(2*pw/3.1416),zero,5e-6,rof1);
        power(v10,TODEV);
       }
      else
       {
        if ((d2/2)>0)
         delay(d2/2 -1.0e-6 -rof1 -(2*pw/3.1416));
        power(v13,TODEV);
        shaped_pulse(selshape,selpw,zero,rof1,rof1); 
        power(v10,TODEV);
        if ((d2/4)>0)
         delay(d2/4 -2*rof1 -2*pw);
        rgpulse(pw,zero,rof1,0.0);
        rgpulse(2*pw,one,0.0,0.0);
        rgpulse(pw,zero,0.0,rof1);
        if ((d2/4)>0)
         delay(d2/4 -2*rof1 -2*pw -2*pw/3.1416);
       }
      rgpulse(pw, v2, rof1, 0.0);
      rgpulse(pw, v3, 1.0e-6, rof2);
   status(C);
}
