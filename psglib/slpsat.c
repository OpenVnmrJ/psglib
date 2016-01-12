/* slpsat - standard two pulse sequence with optional composite
            observe pulse and SLP presaturation with xmtr

         d1 - satdly - p1 - d2 - pw - at

   satflg = 'y' : obs xmtr saturation at satfrq with power satpwr
          = 'n' : dec saturation at dof with power dpwr
   composit = 'y' : uses composite 90 for pw
            = 'n' : uses normal 90 for pw   
   slpflg = 'y' : slp presaturation 
   shape_pw = length of the individual shaped_pulse element
                (calculated from h2off)
   h2off  = water off-resonance frequency (typically tof - satfrq)
   satshape = shape of the shaped_pulse element (use slpsatd or slpsatu)

---------------------------------------------------------------------------
                satshape = slpsatd (positive h2off)
                           slpsatu (negative h2off)
           slpsatd is a 72 step hard pulse shape with the phase of
           one element differs from the previous one by -50 degrees
           slpsatu is a 72 step hard pulse shape with the phase of
           one element differs from the previous one by +50 degrees

           slpsatu.RF and slpsatd.RF must exist in the shapelib


   HINT :  After setting the h2off and satshape array tof (typically
           +/- 10 Hz from satfrq+h2off) to find the best suppression.
           This is to compensate for the "round-off" errors in the 
           shape_pw calculation.
----------------------------------------------------------------------------

         -   KRISH Dec. 12, 1991 
             (revised)  Jan 13, 1992                                 */


#include <standard.h>

static int phs1[8] = {0,2,1,3,0,2,1,3},
           phs2[8] = {1,1,2,2,3,3,0,0};

pulsesequence()
{
   double satdly,
          satpwr,
          cycles,
          shape_pw,
          h2off,
          satfrq;
   char   composit[MAXSTR],
          slpflg[MAXSTR],
          satshape[MAXSTR],
          satflg[MAXSTR];

   satdly=getval("satdly");
   satfrq=getval("satfrq");
   satpwr=getval("satpwr");
   h2off = getval("h2off");
   getstr("composit",composit);
   getstr("satflg",satflg);
   getstr("satshape",satshape);
   getstr("slpflg",slpflg);

/* CHECK HARDWARE CONFIGURATION */

   if ((rfwg[0] != 'y') && (slpflg[0] == 'y'))
    {
     printf("NOTE: slpflg=y  requires PPM in Channel 1. P.S.G. aborted.\n");
     abort(1);
    }

   if ((slpflg[0]=='y') && (h2off != 0.0))
   {
   shape_pw = 10/h2off;
   if (shape_pw < 0.0)
    shape_pw = -shape_pw;

   cycles = (satdly/(shape_pw + 15.4e-6));
   cycles = 2.0*(double) (int) (cycles/2.0);
 
   initval(cycles,v11);
   }


   settable(t1,8,phs1);
   settable(t2,8,phs2);

   getelem(t1,ct,oph);

   add(oph,one,v4);
   add(oph,two,v5);
   add(oph,three,v6);


   status(A);

	if (satflg[0] == 'y')
          {
           rlpower(satpwr,TODEV);
           hsdelay(d1);
           if (slpflg[0] == 'n')              /* standard presat */
           {
            if (satfrq != tof)
	     offset(satfrq,TODEV);
	    rgpulse(satdly,t2,40.0e-6,0.0);
            if (satfrq != tof)
	     offset(tof,TODEV);
	   }
           if (slpflg[0] == 'y')              /*  slp presat  */
           {
            if (h2off !=0.0)
            {
            rcvroff();
            delay(rof1);
            starthardloop(v11);
             shaped_pulse(satshape,shape_pw,t2,0.0,1e-6);
            endhardloop();
            delay(rof2);
            rcvron();
            }
            else                             /* presat with WFG */
            {
             if (satfrq != tof)
              offset(satfrq,TODEV);
             shaped_pulse(satshape,satdly,t2,40.0e-6,0.0);
             if (satfrq != tof)
              offset(tof,TODEV);
            }
           }
           rlpower(tpwr,TODEV);
           delay(40.0e-6);
          }
	 else
	  hsdelay(d1+satdly);

   status(B);

	pulse(p1,zero);
	hsdelay(d2);


	if (composit[A] == 'n')
	   rgpulse(pw,oph,rof1,rof2);
	else
	 {
	  rgpulse(pw,v4,rof1,1.0e-6);
	  rgpulse(pw,v5,0.0,1.0e-6);
	  rgpulse(pw,v6,0.0,1.0e-6);
	  rgpulse(pw,oph,0.0,rof2);
	 }
   status(C);
}

