
/* wnoesy -    NOESY wherein each pulse may be individually adjusted
               for power and shape. Useful for selective NOESY (allowing
               reduction of sw1 and ni). Use lsfrq1 to move folded peaks
               in F1.

                ....d1....satdly..p1(SLP)..d2...p2..mix..p3..at

                     p1,p2 and p3 use wfg shaped pulses
                     p1 pulse is SLP(shifted laminar pulse) to provide
                     a 90 degree pulse in region of interest. A shapelib.RF
                     file "p1shape" must be created prior to use
                     via the "convolute" macro. The rf carrier is always
                     centered on the water at "tof".

  Parameters:

       mix = mixing time.
     phase =   0: gives non-phase-sensitive experiment (P-type peaks);
                  nt(min) = multiple of 16
                  nt(max) = multiple of 64

             1,2: gives HYPERCOMPLEX phase-sensitive experiment;
               3: gives TPPI phase-sensitive experiment;
                  nt(min) = multiple of  8
                  nt(max) = multiple of 32

    satmode= 'y' uses obs xmtr at power satpwr and tof
              satmode='yyyn', for example to saturate during satdly,d2 and mix
    satdly = length of presaturation period;
    p1pwr  = power level for first(tophat) pulse using shape "p1shape"
    p2pwr  = power level for second pulse using shape "p2shape"
    p3pwr  = power level for third pulse using shape "p3shape"
   satpwr  = saturation power level
             G.Gray  Palo Alto  Sept 1991 
*/


#include <standard.h>

pulsesequence()
{
   double          arraydim,
                   ss,
                   satpwr,p1pwr,p2pwr,p3pwr,
                   satdly,
                   ni,p2,p3,
                   mix;
   int             t1_counter,
                   iphase;
   char            p1shape[MAXSTR],p2shape[MAXSTR],p3shape[MAXSTR],
                   satmode[MAXSTR];


/* LOAD VARIABLES */
   ni = getval("ni");
   arraydim = getval("arraydim");
   satdly = getval("satdly");
   p2 = getval("p2"); p3 = getval("p3");
   satpwr = getval("satpwr");
   p1pwr  = getval("p1pwr");
   p2pwr  = getval("p2pwr");
   p3pwr  = getval("p3pwr");
   mix = getval("mix");
   iphase = (int) (getval("phase") + 0.5);
   ss = getval("ss");
   getstr("satmode", satmode);
   getstr("p1shape", p1shape);
   getstr("p2shape", p2shape);
   getstr("p3shape", p3shape);
   initval(p1pwr,v13);
   initval(p2pwr,v4);
   initval(satpwr,v8);
   initval(p3pwr,v11); 

   if (iphase == 3)
   {
      t1_counter = ((int) (ix - 1)) / (arraydim / ni);
      initval((double) (t1_counter), v14);
   }
   else
      assign(zero, v14);


/* CHECK CONDITIONS */
   if ((rof1 < 9.9e-6) && (ix == 1))
      fprintf(stdout,"Warning:  ROF1 is less than 10 us\n");

   if (satpwr > 40)
        {
         printf("satpwr too large - acquisition aborted./n");
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
   if (iphase == 0)
   {
      assign(v10, v9);
      hlv(v10, v10);
      mod2(v9, v9);
   }
   else
   {
      assign(zero, v9);
   }
   assign(v10,v1);
   hlv(v10, v10);
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
   add(v10,v14,v5);
   if (iphase == 2)
      { incr(v2); incr(v5); }
   if (iphase == 3)
      add(v2, v14, v2);		/* TPPI phase increment */


/*HYPERCOMPLEX MODE USES REDFIELD TRICK TO MOVE AXIAL PEAKS TO EDGE */
    initval(2.0*(double)(((int)(d2*getval("sw1")+0.5)%2)),v6);
  if ((iphase==1)||(iphase==2))
       {add(v2,v6,v2); add(oph,v6,oph); add(v5,v6,v5);}  

/* BEGIN THE ACTUAL PULSE SEQUENCE */
   status(A);
      power(zero,DODEV);
      hsdelay(d1);
      if (satmode[A] == 'y')
      {
        power(v8,TODEV);
        rgpulse(satdly,v5,rof1,rof1);
      }
   status(B);
      power(v13,TODEV);
      shaped_pulse(p1shape,p1, v2, rof1, 1.0e-6); /* use SLP pulse */
      if (satmode[B] =='y')
       {
        if (d2 > 0.0)
         {
          power(v8,TODEV);
          rgpulse(d2-15.6e-6-rof1,zero,5.0e-6,5.0e-6);
         }
       }
      else
       {  
        if (d2 > 0.0)
         delay(d2-11.0e-6);
       }
      power(v4,TODEV);
      shaped_pulse(p2shape,p2, v1, 10.0e-6, 1.0e-6);
   status(C);
      if (satmode[C] == 'y')
       {
         power(v8,TODEV);
          if (hs[C] == 'y')
           {
             hsdelay(hst);
             rgpulse(mix-hst,zero,20.0e-6,rof1);
           }
          else rgpulse(mix,zero,20.0e-6,rof1);
       }
        else hsdelay(mix);
      
   status(D);
      power(v11,TODEV);
      shaped_pulse(p3shape,p3, v3, rof1, rof2);
}
