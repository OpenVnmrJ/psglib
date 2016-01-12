
/* goesy.c            DPFGSENOE with selective excitation 

                    bbpwr    = sufficient for full inversion of 1H spectrum
                    tpwr     = normal 90 power
                    bbpw     = pw for sech180  (should cover 12 ppm)
                    selpwr   = selective 180 power
                    selshape  (selective shape)
                    bb180    = name of broadband inversion pulse (typ.sech180) 
                  Some typical values:
                    gzlvl1   = -500
                    gzlvl2   = -3000
                    gzlvl3   = +500
                    gtime1   = 1 ms
                    gtime2   = 5 ms
                    gtime3   = 2 ms
                    mix      = desired mixing time, say 500-1000 ms
                               for small molecule
                    gstab1   = 500 us or so
                    gstab2   = 500 us or so

Set TOF to the center of the spectral window. 
use Pbox to select target multiplet and shape
*/

#include <standard.h>

pulsesequence()
{

   double         
                   mix,
                   gzlvl1,
                   gzlvl2,
                   gzlvl3,
                   gtime1,
                   gtime2,
                   gtime3,
                   gstab2,
                   gstab1,
                   selpwr,
                   selpw,
                   bbpw,
                   bbpwr,satpwr;
   char            selshape[MAXSTR];    
   char            bb180[MAXSTR];   /*  bb180 is sech180.RF */

   satpwr = getval("satpwr");
   selpw = getval("selpw");
   bbpwr = getval("bbpwr");
   bbpw = getval("bbpw");
   selpwr  = getval("selpwr");
   mix = getval("mix");
   gzlvl1 = getval("gzlvl1");
   gzlvl2 = getval("gzlvl2");
   gzlvl3 = getval("gzlvl3");
   gtime1 = getval("gtime1");
   gtime2 = getval("gtime2");
   gtime3 = getval("gtime3");
   gstab1 = getval("gstab1");
   gstab2 = getval("gstab2");

    getstr("selshape",selshape);
    if (selshape[0] == '\0')
    {
      text_error("no selshape? ABORT");
      abort(1);
    }
   getstr("bb180",bb180);

/* STEADY-STATE PHASECYCLING  */


/* CALCULATE PHASES */

     mod4(ct,v2);
     hlv(ct,v1);
     hlv(v1,v1);
     mod4(v1,v1);
     dbl(v2,oph);
     dbl(v1,v9);
     add(v9,oph,oph);

/* CYCLOPS */

     hlv(ct,v10);
     hlv(v10,v10);
     hlv(v10,v10);
     hlv(v10,v10);
     mod4(v10,v10);
     add(one,v10,v4);

     add(v10,v1,v1);
     add(v10,v2,v2);
     add(v10,oph,oph);

   /* equilibrium period */

status(A);
     delay(d1);
     obspwrf(4095.0);
     obspower(satpwr);
     rgpulse(satdly,zero,rof1,rof2);
     obspower(tpwr);

   /* hard 90, G1-selective 180-G1,G2-selctive 180-G2 (to refocus the
target)  */

status(B);

     rgpulse(pw,v10,rof1,rof2);

     obspower(selpwr);
     delay(gstab1);
     zgradpulse(gzlvl1,gtime1);
     delay(gstab1);

     shaped_pulse(selshape,selpw,v2,rof1,rof2);

     delay(gstab1);
     zgradpulse(gzlvl1,gtime1);
     delay(gstab1);

     delay(gstab1);
     zgradpulse(gzlvl2,gtime1);
     delay(gstab1);

     shaped_pulse(selshape,selpw,v1,rof1,rof2);

     delay(gstab1);
     zgradpulse(gzlvl2,gtime1);
     delay(gstab1);
     obspower(tpwr);
     obspwrf(4095.0);

 status(C);

     rgpulse(pw,v10,rof1,rof2);
     zgradpulse(gzlvl3,gtime2);
     obspower(bbpwr);
  /* mixing along the z axis, then dephase any signal at xy plane */
   if ( mix > 0.0)
    {
       delay((31.0*mix/100.0)-gtime3-gstab2-gtime2-(bbpw/2.0));
       zgradpulse(gzlvl3,gtime3);
       delay(gstab2);
       shaped_pulse(bb180,bbpw,zero,rof1,rof2);
       delay(gstab2);
       zgradpulse(-gzlvl3,gtime3);
       delay((49.0*mix/100.0)-gtime3*3-gstab2*2-bbpw);
       zgradpulse(gzlvl3,gtime3*2.0);
       delay(gstab2);
       shaped_pulse(bb180,bbpw,zero,rof1,rof2);
       delay(gstab2);
       zgradpulse(-gzlvl3,gtime3*2.0);
       delay((20.0*mix/100.0)-gtime3*2.0-gstab2-(bbpw/2.0));
    }
status(D);
     obspower(tpwr);
     rgpulse(pw,v10,rof1,rof2);
}



