/* hcosyinv - heteronuclear COSY with proton observe 


   H1  -      BB  -               - 90 - Acq.
   X          d1  - 90 -    t1    - 90

   pwxlvl - X pulse power level
   pwx    - X 90 deg pulse
   hdpwr  - proton BB decoupling (during d1 using obs channel) power 
   hdres  - tip-angle resolution for the hdshape
   hd90   - 90 pulse width at hdpwr for BB decoupling
   hdshape - decoupling modulation shape
   hdflg - 'y' for decoupling
           'n' no decoupling
   phase - 1,2 (only hypercomplex with FAD is supported)
   sspul - does trim(x)-trim(y) before d1 */

#include <standard.h>

pulsesequence()
{
  double   pwxlvl,
           pwx,
           hdpwr,
           hdres,
           hd90,
           phase;

  char     hdshape[MAXSTR],
           sspul[MAXSTR],
           hdflg[MAXSTR];

  int      iphase;


  pwxlvl = getval("pwxlvl");
  pwx    = getval("pwx");
  hdpwr  = getval("hdpwr");
  hd90   = getval("hd90");
  hdres  = getval("hdres");
  phase  = getval("phase");

  getstr("hdshape",hdshape);
  getstr("hdflg",hdflg);
  getstr("sspul",sspul);

  iphase = (int)(phase + 0.5);

  hlv(ct,v1);                    /* v1 = 0 0 1 1 2 2 3 3 */
  mod4(ct,v1);

  mod2(ct,v2);                  /*  v2 = 0 1 0 1 0 1 0 1 */
  dbl(v2,v2);                   /*  v2 = 0 2 0 2 0 2 0 2 */


  add(v2,v1,v2);                /*  v2 = 0 2 1 3 2 0 3 1 */

  assign(v2,oph);

  initval(2.0*(double)(((int)(d2*getval("sw1")+0.5)%2)),v14);

  if (iphase == 2)
   incr(v2);

  add(v2,v14,v2);
  add(oph,v14,oph);

  status(A);

   rlpower(pwxlvl,DODEV);
   rlpower(tpwr,TODEV);

   if (sspul[0] == 'y')
    { hsdelay(hst);
      rgpulse(pw,v1,rof1,rof2);
      hsdelay(hst+.001);
    }

   if (hdflg[0] == 'y')
    { rlpower(hdpwr,TODEV);
      delay(0.1);
      rcvroff();
      obsprgon(hdshape,hd90,hdres);
      xmtron();
      delay(d1);
      xmtroff();
      obsprgoff();
      rlpower(tpwr,TODEV);
      delay(0.05);
      rcvron();
    }
   else
    delay(d1);
    rcvroff();
    decphase(v2);

  status(B);
    decrgpulse(pwx,v2,rof1,1e-6);
    if (d2 > 0.0)
     delay(d2 - 2e-6 - (4*pwx/3.1416));
    else
     delay(d2);
    simpulse(pw,pwx,v1,v1,1e-6,rof2);

  status(C);
    rcvron();
}
