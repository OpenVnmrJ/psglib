/* hmqcnoesy3rf.c - 

        supports both hypercomplex (Bax mode with shifted axial peaks) 
	and TPPI implementations, steady-state
        transients for either all or just the first t1 increment,
        pulse+receiver phasecycling during steady-state transients, and
        presaturation using the transmitter.

        USES THIRD CHANNEL FOR X PULSES.  1H IN F2, X IN F1 CORRELATION

   Parameters:

     satmode = 'yn':  presaturation during relaxation period (satdly) with xmtr
              'yyy':  presaturation during relaxation,(optional null) and mix time
     satfrq = presaturation frequency
     satdly = saturation time during the relaxation period
     satpwr = saturation power for all periods of presaturation with xmtr
        dof = frequency for decoupling during acquisition 
         hs = 'yn':  homospoil pulse (hst) during the d1 relaxation delay
              'yy':  homospoil pulse during both d1 and nulling period (null)
         pw = 90 degree transmitter pulse length for protons (power = TPWR)
     pwx2lvl = dpwr2 level for decoupler pulses
        pwx2 = 90 degree decoupler pulse length for the heteronucleus
          j = one-bond heteronuclear coupling constant (in Hz)
         dm2 = 'nny':  activates 1/(2j) refocusing period before acquisition
              and includes heteronuclear broadband decoupling
                       (recommended)
      phase = 1,2: hypercomplex phase-sensitive experiment with F1 quadrature
			(using axial shifting phase cycling of Bax);
                3: TPPI phase-sensitive experiment with F1 quadrature
       null = nulling time for protons not attached to observed heteronucleus
         nt = multiple of 8
	mix = mixing time 

     Contact-  G. Gray  (palo alto)
     revision-   modified from hmqc.c; modified 1/13/93 to fix d2 correction
                 to use pwx2 instead of pwx
Program aborts if:
	at >= 250 msec
	Duty cycle > 20%
	dm2 = 'yyy','yny','nyy','yyn',. . .*/

#include <standard.h>

pulsesequence()
{
/* VARIABLE DECLARATION */
    double ss,at,pwx2lvl,phase,pwx2,j,bird,null,satfrq,satdly,satpwr;
    double mix,ni,arraydim;
    int iphase;
    char satmode[MAXSTR];


/* LOAD VARIABLES */
    ni = getval("ni");
    arraydim = getval("arraydim");
    satfrq = getval("satfrq");
    satdly = getval("satdly");
    satpwr = getval("satpwr");
    pwx2lvl = getval("pwx2lvl");
    mix = getval("mix");
    at = getval("at");
    phase = getval("phase");
    pwx2 = getval("pwx2");
    null = getval("null");
    j = getval("j");
    ss = getval("ss");
    getstr("satmode",satmode);


/* CHECK CONDITIONS */

    if (at >= .250)
	{
	printf("at >= 250msec.  Please reduce. \n");
	abort(1);
	}
    if ((dm2[A] == 'y') || (dm2[A] == 'y') )
	{
	printf("dm2 must be nnn or nny. \n");
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
	
    initval(pwx2lvl,v10); initval(dpwr2,v8);
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
   initval(2.0*(double)((int)(d2*getval("sw1")+0.5)%2),v5);
   if ((iphase==1)||(iphase==2)) 
        {add(v3,v5,v3); add(oph,v5,oph);}
    if (iphase == 3) 
	 add(v9,v3,v3);             /* TPPI phase increment */


/* BEGIN ACTUAL PULSE SEQUENCE CODE */
    status(A);
       power(v10,DO2DEV);
       power(v6,TODEV);
       hsdelay(d1);

/* selective saturation period */
       if (satmode[A] == 'y')
         {
          offset(satfrq,TODEV);
          power(v7,TODEV);
          rgpulse(satdly,zero,4.0e-5,1.0e-5);
          offset(tof,TODEV);
          power(v6,TODEV); 
         }

/* if null is 0 eliminate bird inversion pulse */
/* bird pulse - 180 for protons not coupled to 13C */

    status(B);
       if (null != 0.0)
         {
          rgpulse(pw,v1,rof1,rof2); 
          delay(bird);
          dec2rgpulse(pwx2,v4,rof1,rof1);
          simpulse(2.0*pw,2.0*pwx2,v1,v1,rof1,rof1); 
          dec2rgpulse(pwx2,v4,rof1,rof1);
          delay(bird);
          rgpulse(pw,v2,rof1,rof2); 

/* nulling time for protons not coupled to 13C */

          if (satmode[B] == 'y')
            {
             offset(satfrq,TODEV);
             power(v7,TODEV);
             rgpulse(null,zero,4.0e-5,1.0e-5);
             offset(tof,TODEV);
             power(v6,TODEV); 
            }
          else hsdelay(null);
         }

       add(one,v13,v2);
       rcvroff();
       rgpulse(pw,v1,rof1,rof2);
      

           delay(bird);
	   dec2rgpulse(pwx2,v3,rof1,rof1);
           if (d2/2.0 > 0.0)
           delay(d2/2.0 - 2*rof1 - 2*pwx2/3.1414 - pw);
           else
           delay(d2/2.0);
           rgpulse(2.0*pw,v1,rof1,rof1);
           if (d2/2.0 > 0.0)
           delay(d2/2.0 - 2*rof1 - 2*pwx2/3.1414 - pw);
           else
           delay(d2/2.0);
	   dec2rgpulse(pwx2,v1,rof1,rof1);

       power(v8,DO2DEV);
       dec2phase(zero);          /* to allow x decoupling */
       if (dm2[C] == 'y')
	 delay(bird);
       rgpulse(pw,v1,rof1,0.0); 
       if (satmode[B] == 'y')
         {
          offset(satfrq,TODEV);
          power(v7,TODEV);
          rgpulse(mix,zero,4.0e-5,1.0e-5);
          offset(tof,TODEV);
          power(v6,TODEV);
         }
       else
         delay(mix);
    status(C);
       rgpulse(pw,v1,rof1,rof2);
       txphase(zero); rcvron();
}
