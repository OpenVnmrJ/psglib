/* hmqctocsy- HMQC with CLEAN-TOCSY(MLEV17) proton relay;
            supports both hypercomplex and TPPI implementations, steady-state
            transients for either all or just the first t1 increment,
            pulse+receiver phasecycling during steady-state transients, and
            presaturation using the transmitter.


   Parameters:

     window = clean-tocsy window(in us)
     satmode = 'yn':  presaturation during relaxation period (satdly) with xmtr
              'yy':  presaturation during both relaxation period (satdly)
                     and null period (null) with xmtr
     satfrq = presaturation frequency
     satdly = saturation time during the relaxation period
     satpwr = saturation power for all periods of presaturation with xmtr
        dof = frequency for decoupling during acquisition 
         hs = 'yn':  homospoil pulse (hst) during the d1 relaxation delay
              'yy':  homospoil pulse during both d1 and the nulling period (null)
         p1 = 90 degree transmitter pulse length for protons (power = p1lvl)
       tpwr = power level for MLEV-16 spin lock and trim pulses 
         pw = 90 degree transmitter pulse for the MLEV-16 spin lock
     pwxlvl = dpwr level for decoupler pulses
        pwx = 90 degree decoupler pulse length for the heteronucleus
          j = one-bond heteronuclear coupling constant (in Hz)
        mix = mixing time for isotropic mixing
       trim = length of trim pulse preceding the MLEV-16 spin lock
         dm = 'nny':  activates 1/(2j) refocusing period before acquisition
                for heteronuclear broadband decoupling
      phase = 1,2: hypercomplex phase-sensitive experiment with F1 quadrature
                3: TPPI phase-sensitive experiment with F1 quadrature
       null = nulling time for protons not attached to observed heteronucleus
         nt = multiple of 8


   s. farmer    25  february    1988
   g. gray      14  february    1990 */
#include <standard.h>
#define PSG_LC
#include <lc.h>
mleva()
{
 rgpulse(pw,v11,0.0,0.0); delay(getval("window"));
 rgpulse(2.0*pw,v12,0.0,0.0); delay(getval("window"));
 rgpulse(pw,v11,0.0,0.0);
}

mlevb()
{
 rgpulse(pw,v13,0.0,0.0); delay(getval("window"));
 rgpulse(2.0*pw,v2,0.0,0.0); delay(getval("window"));
 rgpulse(pw,v13,0.0,0.0);
}

pulsesequence()
{
/* VARIABLE DECLARATION */
    double ss,pwxlvl,phase,pwx,j,bird,null,satfrq,satdly,satpwr;
    double cycles,trim,p1lvl,mix,window;
    char tn[MAXSTR],satmode[MAXSTR];
    int iphase;


/* LOAD VARIABLES */
    window = getval("window");
    satfrq = getval("satfrq");
    satdly = getval("satdly");
    satpwr = getval("satpwr");
    pwxlvl = getval("pwxlvl");
    mix = getval("mix");
       ss = getval("ss");
    p1lvl = getval("p1lvl");
    trim = getval("trim");
    dpwr = getval("dpwr");
    phase = getval("phase");
    pwx = getval("pwx");
    null = getval("null");
    j = getval("j");
    ss = getval("ss");
    getstr("tn",tn);
    getstr("satmode",satmode);


/* INITIALIZE VARIABLES */
    if (j > 0.0) bird = 1.0/(2.0*j); else bird = 0.0;
    iphase = (int) (phase + 0.5);
    if (iphase == 3) initval((double) (ix-1),v14);
    initval(pwxlvl,v10); initval(dpwr,v8); initval(p1lvl,v9);
    initval(satpwr,v7); initval(tpwr,v6);
    if (ss < 0)
      {
       ss = (-1)*ss; initval(ss,ssctr);
      }
    else if ((ss > 0) && (ix == 1)) ss = ss;
    initval(ss,ssval);


/* PHASECYCLE CALCULATION */
    ifzero(ssctr);
       hlv(ct,v1); mod2(ct,v3);
    elsenz(ssctr);
       sub(ssval,ssctr,v13);
       hlv(v13,v1); mod2(v13,v3);
    endif(ssctr);

    hlv(v1,v11); hlv(v11,v11);
    add(v1,one,v4);                             /* for composite 180 */
    add(two,v1,v2); 
    dbl(v3,v3); add(v3,v1,v3); assign(v3,oph);
    mod2(v11,v11); dbl(v11,v11); add(v11,v1,v11);
    add(one,v11,v12); add(one,v12,v13);
    if (iphase == 2) incr(v3);
    if (iphase == 3) add(v14,v3,v3);             /* TPPI phase increment */


/* BEGIN ACTUAL PULSE SEQUENCE CODE */
    status(A);
       delay(1.0e-5); power(v10,DODEV);
       power(v9,TODEV);
       rcvroff(); hsdelay(d1);

/* selective saturation period */
       if (satmode[A] == 'y')
         {
          offset(satfrq,TODEV);
          power(v7,TODEV);
          rgpulse(satdly,zero,4.0e-5,1.0e-5);
          offset(tof,TODEV);
          power(v9,TODEV); rcvroff();
          delay(5.0e-5);
         }

/* if null is 0 eliminate bird inversion pulse */
/* bird pulse - 180 for protons not coupled to 13C */

    status(B);
       if (null != 0.0)
         {
          rgpulse(p1,v1,1.0e-5,0.0); rcvroff();
          delay(bird);
          decpulse(pwx,v4);
          simpulse(2.0*p1,2.0*pwx,v1,v1,0.0,0.0); rcvroff();
          decpulse(pwx,v4);
          delay(bird);
          rgpulse(p1,v2,0.0,0.0); rcvroff();

/* nulling time for protons not coupled to 13C */

          if (satmode[B] == 'y')
            {
             offset(satfrq,TODEV);
             power(v7,TODEV);
             rgpulse(null,zero,4.0e-5,1.0e-5);
             offset(tof,TODEV);
             power(v9,TODEV); rcvroff();
             delay(5.0e-5);
            }
          else hsdelay(null);
         }


       add(one,v13,v2);
       rgpulse(p1,v1,0.0,0.0); rcvroff();
       delay(bird); decpulse(pwx,v3);
       delay(d2/2.0);
       rgpulse(2.0*p1,v1,0.0,0.0); rcvroff();
       delay(d2/2.0); decpulse(pwx,v1);
       if (dm[C] == 'y') delay(bird);
       power(v8,DODEV); txphase(v12);
       power(v6,TODEV);
       decphase(zero);          /* to allow x decoupling */
  status(C);
/* calculate and initialize loop counter */
       cycles = (mix-trim)/(64.66*pw+32*window);
       cycles = 2.0*(double)(int)(cycles/2.0);
       initval(cycles,v5);
       rgpulse(trim,v12,5.0e-6,0.0);
       if (cycles > 1.0)
         {
          starthardloop(v5);
            mleva(); mlevb(); mlevb(); mleva();
            mlevb(); mlevb(); mleva(); mleva();
            mlevb(); mleva(); mleva(); mlevb();
            mleva(); mleva(); mlevb(); mlevb();
            rgpulse(0.66*pw,v12,0.0,0.0);
          endhardloop();
         }
       txphase(zero); delay(rof2); rcvron();
  status(D);
}
