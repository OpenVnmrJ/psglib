/* w1dtocsyz   -  1D tocsy for rapid scalar coupling correlations;
                   internal subtraction of on-resonance and off-
                   resonance selective excitation.
                   Uses wfg-shaped pulse inversion of desired resonance on 
                   alternate scans, followed by normal tocsy sequence.
                   Uses wfg-shaped solvent presaturation
                   Uses "clean-tocsy" via window parameter
                   Uses optional z-filtering
     USES OFFSET TO MOVE TRANSMITTER FREQUENCY TO SATFRQ,SELFRQS,AND TOF

    ref: a. bax and d.g. davis, j. magn. reson. 65, 355 (1985)

    pw - 90 degree pulse during mlev periods (at power level tpwr)
    p1 - 90 degree excitation pulse (at power p1pwr)
    sspul- 'y' kills all magnetization by trim(x)-trim(y)
    sattime - presaturation delay using obs. xmtr.
    satfrq - presaturation frequency
    satpwr - presaturation power level
    satshape - type of presat pulse shape (must be in shapelib)
    intsub - 'y': internal subtraction of data acquired by on-resonance
                  and off-resonance selective excitation
             'n': data acquired by on-resonance and off-resonance selec-
                  tive excitation are stored separately
    selfrqon -  frequency of transmitter selective pulse (on resonance);
                array if intsub = 'n'
    selfrqof - frequency of transmitter selective pulse (off resonance);
                an inactive parameter if intsub = 'n'
    selpwr - power of transmitter selective pulse
             (remember to add 20db to selpwr if using a freq.-shifted pulse)
    selshape  -  shape of selective inversion pulse (must be in shapelib)
    seltime - length(in ms) of transmitter selective (180) pulse at frequency satfrq.
               (create seltime parameter as a "pulse")
    mix - mixing time
    window - "clean tocsy" window parameter(in us)  (create parameter as pulse)
    nt - multiple of 32  (intsub = 'n')
         multiple of 64  (intsub = 'y')

                  G.Gray  Palo Alto   Sept 1991 */


#include <standard.h>

mleva()
{
 double window;

 window = getval("window");
 txphase(v2); delay(pw);
 xmtroff(); delay(window); xmtron();
 txphase(v3); delay(2.0*pw);
 xmtroff(); delay(window); xmtron();
 txphase(v2); delay(pw);
}

mlevb()
{
 double window;

 window = getval("window");
 txphase(v4); delay(pw);
 xmtroff(); delay(window); xmtron();
 txphase(v5); delay(2.0*pw);
 xmtroff(); delay(window); xmtron();
 txphase(v4); delay(pw);
}


pulsesequence()
{
/*DEFINE LOCAL VARIABLES */
  double satfrq,satpwr,sattime,p1pwr,mix,cycles,window,selfrqon,
         selfrqof,selpwr,seltime,zfilter;
    char sspul[MAXSTR],intsub[MAXSTR],satshape[MAXSTR],selshape[MAXSTR];


/* LOAD AND INITIALIZE VARIABLES */
  selfrqon = getval("selfrqon");
  selfrqof = getval("selfrqof");
  getstr("selshape",selshape);
  seltime = getval("seltime");
  selpwr = getval("selpwr");
  mix = getval("mix");
  p1pwr = getval("p1pwr");
  getstr("satshape",satshape);
  getstr("sspul",sspul);
  sattime = getval("sattime");
  satfrq = getval("satfrq");
  satpwr = getval("satpwr");
  zfilter= getval("zfilter");
  window = getval("window"); getstr("intsub",intsub);

  initval(p1pwr,v8);
  initval(tpwr,v7);
  initval(satpwr,v11);
  initval(selpwr,v6);

/* CALCULATE PHASES */
 if (intsub[A] == 'y')
 {
    hlv(ct,v1);
 }
 else
 {
    assign(ct,v1);
 }
 hlv(v1,v13);
 hlv(v13,v13);
 hlv(v13,v2);
 mod2(v13,v13);
 mod4(v1,v1);
 add(one,v1,v1);
 add(v1,v13,v1);
 assign(v1,oph);
 mod4(v2,v2);
 add(one,v2,v2); 
 add(v13,v2,v2);
 add(one,v2,v3);
 add(one,v3,v4);
 add(one,v4,v5);
 if (intsub[A] == 'y')
   {
    mod2(ct,v14);    /* trigger for the alteration of the saturation freq */
    ifzero(v14);
       add(oph,two,oph);
    endif(v14);
   }
/* calculate and initialize loop counter */
       cycles = (mix)/(66.0*pw+32*window);
       cycles = 2.0*(double)(int)(cycles/2.0);
       mix = cycles*(64.66*pw+32*window);
       initval(cycles,v9);

/* CHECK CONDITIONS */
   if ((selpwr>40)||(seltime>1000))
     {printf("selpwr and/or seltime too large.\n"); abort(1);}

/* BEGIN ACTUAL PULSE SEQUENCE CODE */

 status(A);
    power(zero,DODEV);
    if (sspul[A] == 'y')
     { rgpulse(200*p1,zero,rof1,0.0); rgpulse(200*p1,one,0.0,rof2); }
    hsdelay(d1);
 status(B);
     {
      offset(satfrq,TODEV);
      power(v11,TODEV);
      shaped_pulse(satshape,sattime,zero,rof1,rof2); 
      power(v8,TODEV); 
     }

/* selective pulse or transmitter saturation */
    if (seltime>0)
    {
      power(v6,TODEV);       /* sets obs atten floor for shaped pulse  */
      if (intsub[A] == 'n')
      {
         offset(selfrqon,TODEV);
      }
      else
      {
         ifzero(v14);
            offset(selfrqon,TODEV);
         elsenz(v14);
            offset(selfrqof,TODEV);
         endif(v14);
      }
     shaped_pulse(selshape,seltime*1000,zero,rof1,rof1);
     offset(tof,TODEV);
    }
     power(v8,TODEV); 
     rgpulse(p1,v1,30.0e-6,0.2e-6);
     power(v7,TODEV); rcvroff(); delay(8.0e-6);

/* tocsy mixing time */
    if (mix > 0.0)
    {
       txphase(v5); 
       xmtron(); 
       {
          starthardloop(v9);
            mleva(); mlevb(); mlevb(); mleva();
            mlevb(); mlevb(); mleva(); mleva();
            mlevb(); mleva(); mleva(); mlevb();
            mleva(); mleva(); mlevb(); mlevb();
            rgpulse(.66*pw,v3,rof1,0.0);
          endhardloop();
       }
      if (zfilter > 0.0)
      {
       power(v8,TODEV);
       rgpulse(p1,v1,0.0,rof2);
       delay(zfilter);
       rgpulse(p1,v1,rof1,rof2);
      }
     }
 status(C);
}
