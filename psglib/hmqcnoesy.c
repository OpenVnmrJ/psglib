/* hmqcnoesy.c - 
        supports both hypercomplex (Bax mode with shifted axial peaks) 
	and TPPI implementations, steady-state
        transients for either all or just the first t1 increment,
        pulse+receiver phasecycling during steady-state transients, and
        presaturation using the transmitter.

   Parameters:

     satmode = 'ynn':  presaturation during relaxation period (satdly) with xmtr
              'yyn':  presaturation during both relaxation period (satdly)
                     and null period (null) with xmtr
              'yyy':  presaturation during relaxation,null and mix time
     satfrq = presaturation frequency
     satdly = saturation time during the relaxation period
     satpwr = saturation power for all periods of presaturation with xmtr
        dof = frequency for decoupling during acquisition 
         hs = 'yn':  homospoil pulse (hst) during the d1 relaxation delay
              'yy':  homospoil pulse during both d1 and nulling period (null)
         pw = 90 degree transmitter pulse length for protons (power = TPWR)
     pwxlvl = dpwr level for decoupler pulses
        pwx = 90 degree decoupler pulse length for the heteronucleus
          j = one-bond heteronuclear coupling constant (in Hz)
         dm = 'nny':  activates 1/(2j) refocusing period before acquisition
              and includes heteronuclear broadband decoupling
                       (recommended)
      phase = 1,2: hypercomplex phase-sensitive experiment with F1 quadrature
			(using axial shifting phase cycling of Bax);
                3: TPPI phase-sensitive experiment with F1 quadrature
       null = nulling time for protons not attached to observed heteronucleus
         nt = multiple of 8
	mix = mixing time 

     Contact-  G. Gray  (palo alto)
     revision-   modified from hmqc.c
Program aborts if:
	at >= 250 msec
	Duty cycle > 20%
	dm = 'yyy','yny','nyy','yyn',. . .*/

#include <standard.h>

pulsesequence()
{
/* VARIABLE DECLARATION */
    double ss,at,pwxlvl,phase,pwx,j,bird,null,satfrq,satdly,satpwr;
    double mix,ni,arraydim;
    int iphase;
    char satmode[MAXSTR];


/* LOAD VARIABLES */
    ni = getval("ni");
    arraydim = getval("arraydim");
    satfrq = getval("satfrq");
    satdly = getval("satdly");
    satpwr = getval("satpwr");
    pwxlvl = getval("pwxlvl");
    mix = getval("mix");
    at = getval("at");
    phase = getval("phase");
    pwx = getval("pwx");
    null = getval("null");
    j = getval("j");
    ss = getval("ss");
    getstr("satmode",satmode);


/* CHECK CONDITIONS */

    if ((at+mix) >= .250)
	{
	printf("(at+mix) >= 250msec.  Please reduce. \n");
	abort(1);
	}
    if ((dm[A] == 'y') || (dm[A] == 'y') )
	{
	printf("dm must be nnn or nny. \n");
	abort(1);
	}

/* INITIALIZE VARIABLES */
    if (j > 0.0) 
	     bird = 1.0/(2.0*j);
	else bird = 0.0;
    iphase = (int) (phase + 0.5);
/*    if (iphase == 3) 
      initval((double) (ix-1),v9);  */
	if (iphase == 3)
	  initval((double)((int)((ix-1)/(arraydim/ni)+1e-6)),v9);
	else
	assign(zero,v9);
	
    initval(pwxlvl,v10); initval(dpwr,v8);
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
    if (iphase == 2)
	 incr(v3);
   initval(2.0*(double)((int)(d2*getval("sw1"))%2),v5);
   if ((iphase==1)||(iphase==2)) 
        {add(v3,v5,v3); add(oph,v5,oph);}
    if (iphase == 3) 
	 add(v9,v3,v3);             /* TPPI phase increment */


/* BEGIN ACTUAL PULSE SEQUENCE CODE */
    status(A);
       delay(1.0e-5); power(v10,DODEV);
       power(v6,TODEV);
       rcvroff(); hsdelay(d1);

/* selective saturation period */
       if (satmode[A] == 'y')
         {
          offset(satfrq,TODEV);
          power(v7,TODEV);
          rgpulse(satdly,zero,4.0e-5,1.0e-5);
          offset(tof,TODEV);
          power(v6,TODEV); rcvroff();
          delay(5.0e-5);
         }

/* if null is 0 eliminate bird inversion pulse */
/* bird pulse - 180 for protons not coupled to 13C */

    status(B);
       if (null != 0.0)
         {
          rgpulse(pw,v1,1.0e-5,0.0); rcvroff();
          delay(bird);
          decpulse(pwx,v4);
          simpulse(2.0*pw,2.0*pwx,v1,v1,0.0,0.0); rcvroff();
          decpulse(pwx,v4);
          delay(bird);
          rgpulse(pw,v2,0.0,0.0); rcvroff();

/* nulling time for protons not coupled to 13C */

          if (satmode[B] == 'y')
            {
             offset(satfrq,TODEV);
             power(v7,TODEV);
             rgpulse(null,zero,4.0e-5,1.0e-5);
             offset(tof,TODEV);
             power(v6,TODEV); rcvroff();
             delay(5.0e-5);
            }
          else hsdelay(null);
         }

       add(one,v13,v2);
       rgpulse(pw,v1,0.0,0.0);
       rcvroff();

           delay(bird);
	   decpulse(pwx,v3);
           delay(d2/2.0);
           rgpulse(2.0*pw,v1,0.0,0.0);
           delay(d2/2.0);
	   decpulse(pwx,v1);

/*       delay(0.2e-6);    */
       power(v8,DODEV);
       decphase(zero);          /* to allow x decoupling */
/*       delay(4.6e-5);            to account for power 5.0e-5 time */
       if (dm[C] == 'y')
	 delay(bird);
       rgpulse(pw,v1,0.0,0.0); 
       if (satmode[B] == 'y')
         {
          offset(satfrq,TODEV);
          power(v7,TODEV);
          rgpulse(mix,zero,4.0e-5,1.0e-5);
          offset(tof,TODEV);
          power(v6,TODEV); rcvroff();
          delay(5.0e-5);
         }
       else
         delay(mix);
    status(C);
       rgpulse(pw,v1,0.0,0.0);
       txphase(zero); delay(rof2); rcvron();
}
