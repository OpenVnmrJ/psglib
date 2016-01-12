/* slpnoesy -  2D cross relaxation experiment WITH SLP PRESATURATION option. 
            It can be performed in
            either phase-sensitive or absolute value mode.  Either
            TPPI or the hypercomplex method can be used to achieve
            F1 quadrature in a phase-sensitive presentation.
  Added features : 1. FAD-hypercomplex
                    2. correction for d2 > 0
                    3. xmtr presaturation 

  Parameters:

       mix = mixing time.
     phase =   0: gives non-phase-sensitive experiment (P-type peaks);
                  nt(min) = multiple of 16
                  nt(max) = multiple of 64

             1,2: gives HYPERCOMPLEX phase-sensitive experiment;
               3: gives TPPI phase-sensitive experiment;
                  nt(min) = multiple of  8
                  nt(max) = multiple of 32

    satflg  = 'y' uses obs xmtr for presat at satfrq and satpwr
              for satdly seconds
                 set dm='n' and satflg='ynyn', for example
    satfrq  = saturation frequency for xmtr presaturation
    satpwr  = saturation power for xmtr presaturation
    satdly  = saturation period follows D1
    slpflg  = flag that turns on slp presaturation
    satshape = shape used for slp presaturation
    satshape = slpsatd (tof > satfrq)
               slpsatu (tof < satfrq)
        ss = number of steady state pulses; if ss < 0 then -ss
             steady-state transients are performed before EVERY increment;
             if ss > 0, then ss steady-state transients are performed
             before only the first increment.
     sspul = 'y': selects for trim(x)-trim(y) sequence at start of pulse sequence
             'n': normal NOESY experiment


      KRISH    12 Dec       1991
      revised  13 Apr 1992 (vvk)  */


#include <standard.h>

pulsesequence()
{
   double          arraydim,
                   ss,
                   satpwr,
                   satdly,
                   satfrq,
                   cycle1,
                   cycle2,
                   shape_pw,
                   h2off,
                   ni,
                   mix;
   int             phase;
   char            sspul[MAXSTR],
                   slpflg[MAXSTR],
                   satshape[MAXSTR],
                   satflg[MAXSTR];

/* LOAD VARIABLES */
   ni = getval("ni");
   arraydim = getval("arraydim");
   satdly = getval("satdly");
   satfrq = getval("satfrq");
   satpwr = getval("satpwr");
   h2off = getval("h2off");
   mix = getval("mix");
   phase = (int) (getval("phase") + 0.5);
   ss = getval("ss");
   getstr("sspul", sspul);
   getstr("satflg",satflg);
   getstr("slpflg",slpflg);
   getstr("satshape",satshape);

   if ((slpflg[0] == 'y') && (h2off != 0.0))
   {
   shape_pw = 10/h2off;
   if (shape_pw < 0.0)
    shape_pw = -shape_pw;

   cycle1 = (satdly/(shape_pw + 15.4e-6));
   cycle1 = 2.0*(double)(int)(cycle1/2.0);
   initval(cycle1,v7);

   cycle2 = (mix/(shape_pw + 15.4e-6));
   cycle2 = 2.0*(double)(int)(cycle2/2.0);
   initval(cycle2,v13);
   }



/* CHECK CONDITIONS */
   if ((rof1 < 9.9e-6) && (ix == 1))
      fprintf(stdout,"Warning:  ROF1 is less than 10 us\n");

   if (rfwg[0] != 'y')
    {
     if ((slpflg[0] == 'y') || (slpflg[2] == 'y'));
    {
     printf("NOTE: slpflg=y requires PPM in Channel 1. P.S.G. aborted.\n");
     abort(1);
    }
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
   assign(v10, v1);
   hlv(v10,v10);
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

  if (phase == 3)                                         /* TPPI */
     initval((double)(ix-1) / (arraydim/ni),v14);
  else
  {
   if ((phase == 1) || (phase == 2))                      /* FAD */
   initval(2.0*(double)(((int)(d2*getval("sw1")+0.5)%2)),v14);
   else
   initval(0.0,v14);
  }

   add(v10,v14,v5);
   if (phase == 2)                                        /* hypercomplex */
     {
      incr(v2);
      incr(v5);
     }

   add(v2, v14, v2);

   if ((phase == 1) || (phase == 2))
     add(oph,v14,oph);
    

/* BEGIN THE ACTUAL PULSE SEQUENCE */
   status(A);
   if (sspul[0] == 'y')
   {
      rgpulse(200*pw, zero, rof1, 0.0);
      rgpulse(200*pw, one, 0.0,rof1);
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
         starthardloop(v7);
          shaped_pulse(satshape,shape_pw,v5,0.0,1e-6);
         endhardloop();
         delay(rof2);
         rcvron();
         }
         else
         {
         if (satfrq != tof)
          offset(satfrq,TODEV);
         shaped_pulse(satshape,satdly,v5,rof1,rof1);
         if (satfrq != tof)
          offset(tof,TODEV);
         }
        }
        else
        {
        if (satfrq != tof)
         offset(satfrq,TODEV);
        rgpulse(satdly,v5,rof1,rof1);
        if (satfrq != tof)
          offset(tof,TODEV);
        }
        rlpower(tpwr,TODEV);
        delay(40.0e-6);
     }

   status(B);

      rgpulse(pw, v2, rof1, 1.0e-6);
      if (d2>0)
       delay(d2-rof1-1.0e-6-(4*pw/3.1416));  /*corrected evolution time */
      else
       delay(d2);
      rgpulse(pw, v1, rof1, 1.0e-6);

   status(C);

      if (satflg[2] == 'y')
        {
         rlpower(satpwr,TODEV);
         if (slpflg[2] == 'y')
         {
          if (h2off != 0.0)
          {
          rcvroff();
          delay(rof1);
          starthardloop(v13);
           shaped_pulse(satshape,shape_pw,zero,0.0,1e-6);
          endhardloop();
          delay(rof2);
          rcvron();        
          }
          else
          {
          if (satfrq != tof)
           offset(satfrq,TODEV);
          shaped_pulse(satshape,mix,zero,2.0e-6,rof1);
          if (satfrq != tof)
           offset(tof,TODEV);
          }
         }
         else
         {
         if (satfrq != tof)
          offset(satfrq,TODEV);
         rgpulse(mix,zero,2.0e-6,rof1);
         if (satfrq != tof)
          offset(tof,TODEV);
         }
         rlpower(tpwr,TODEV);
         delay(40.0e-6);
        }
       else
        hsdelay(mix - rof1);
   status(D);
      rgpulse(pw, v3, rof1, rof2);
}
