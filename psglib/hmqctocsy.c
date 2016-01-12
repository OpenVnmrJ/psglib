/* hmqctocsy- HMQC with CLEAN-TOCSY(MLEV17) proton relay;
            supports both hypercomplex and TPPI implementations, steady-state
            transients for either all or just the first t1 increment,
            pulse+receiver phasecycling during steady-state transients, and
            presaturation using the transmitter.
	    The direct responses are controlled by the parameter mult:
		mult=0 gives both direct & relayed responses inphase
		mult=1 nulls the direct responses. dm='nnnn'- best for this
		mult=2 inverts the phase of the direct responses

 CONTACT:  R. Crouch Burroughs Wellcome Co. (919) 248-3840.

	This sequence was revised from the G. Gray / P. Keifer version.
	The mult business was added and the function of pw & p1 reversed
	so that pw refers to the pulsesequence outside of the mlev time
	and p1 refers to the mlev mix period. All "power" statements
	were replaced with "rlpower".
	Hypercomplex FAD is now supported - highly recommended !!
	A trimX trimY sspulse pulse now supported as well.

   Parameters:

     window = clean-tocsy window(in us)
     satmode = 'yn':  presaturation during relaxation period (satdly) with xmtr
              'yy':  presaturation during both relaxation period (satdly)
                     and null period (null) with xmtr
     satfrq = presaturation frequency
     satdly = saturation time during the relaxation period
     satpwr = saturation power for all periods of presaturation with xmtr
     sstrim = 1000 will give a 1ms trimX,trimY sspulse before d1.
        dof = frequency for decoupling during acquisition 
         hs = 'yn':  homospoil pulse (hst) during the d1 relaxation delay
              'yy':  homospoil pulse during both d1 and the 'null' period
	 hs = 'nn': recommended if sstrim is used. 
         pw = 90 degree transmitter pulse length for protons (power = tpwr)
      p1lvl = power level for MLEV-16 spin lock and trim pulses 
         p1 = 90 degree transmitter pulse for the MLEV-16 spin lock
     pwxlvl = dpwr level for decoupler pulses
        pwx = 90 degree decoupler pulse length for the heteronucleus
          j = one-bond heteronuclear coupling constant (in Hz)
        mix = mixing time for isotropic mixing
       trim = length of trim pulse preceding the MLEV-16 spin lock
      phase = 1,2: hypercomplex phase-sensitive experiment with F1 quadrature
		 with FAD
                3: TPPI phase-sensitive experiment with F1 quadrature
       null = nulling time for protons not attached to observed heteronucleus
         nt = multiple of 8
       mult = multiplier for final carbon pulse:
		 mult=1.0 nulls direct responses. dm='nnnn' recommended. 
                mult=2.0 inverts phase of direct responses.
		mult=0.0 gives normal hmqctocsy; dm='nnyy' dmm='ccss' will
		initiate 13C decoupling during mixing period. 

	newpar: p1lvl-real(63)
		window-pulse(50)
		mix-delay(0.050)
		trim-delay(0.002)
		mult-real(2.0)
		sstrim-pulse(1000)

	NOTE: Sequence aborts if duty cycle >25% or if (at+mix+trim) >= 320msec
	These tests are bypassed if dm='nnnn' as would be best if the
	direct responses were being eliminated. No limit on F2 res under
	these conditions.

   s. farmer    25  february    1988 
   g. gray      14  february    1990 
   p. keifer    15  february    1990
   r. crouch (BW) 20 March 1992 	  */


#include <standard.h>
mleva()
{
 rgpulse(p1,v11,0.0,0.0); delay(getval("window"));
 rgpulse(2.0*p1,v12,0.0,0.0); delay(getval("window"));
 rgpulse(p1,v11,0.0,0.0);
}

mlevb()
{
 rgpulse(p1,v13,0.0,0.0); delay(getval("window"));
 rgpulse(2.0*p1,v2,0.0,0.0); delay(getval("window"));
 rgpulse(p1,v13,0.0,0.0);
}

pulsesequence()
{
/* VARIABLE DECLARATION */
    double at,ss,arraydim,ni,pwxlvl,phase,pwx,j,bird,null,satfrq,satdly;
    double cycles,trim,sstrim,mult,satpwr,p1lvl,mix,window;
    int iphase;
    char satmode[MAXSTR];


/* LOAD VARIABLES */
    window = getval("window");
    satfrq = getval("satfrq");
    satdly = getval("satdly");
    satpwr = getval("satpwr");
    pwxlvl = getval("pwxlvl");
    mix = getval("mix");
    ss = getval("ss");
    at = getval("at");
    p1lvl = getval("p1lvl");
    trim = getval("trim");
    phase = getval("phase");
    pwx = getval("pwx");
    null = getval("null");
    j = getval("j");
    getstr("satmode",satmode);
    sstrim = getval("sstrim");	/* trim X trim Y ss */
    arraydim = getval("arraydim");
    ni = getval("ni");
    mult = getval("mult");     /* Direct response editting control */
 
/* INITIALIZE VARIABLES */
    if (j > 0.0) 
	bird = 1.0/(2.0*j); 
      else bird = 0.0;
    iphase = (int) (phase + 0.5);

    if (iphase == 3)
	 initval((double) (ix - 1) / (arraydim / ni), v14);

    if (ss < 0)
      {
       ss = (-1)*ss; initval(ss,ssctr);
      }
    else if ((ss > 0) && (ix == 1)) ss = ss;
    initval(ss,ssval);

	/* SAFETY CHECKS */

if ((dm[A]=='y') || (dm[B]=='y'))
	{
         printf("dm should equal 'nnnn' or 'nnny' or 'nnyy'.\n");
         abort(1);
        }

	/* Only check duty cycle if we are 13C decoupling */

if ((dm[C]=='y') || (dm[D]=='y'))
 {	
   if (satmode[A] == 'y')
      {
     if (((at+mix+trim)/(at+mix+trim+d1+satdly)) >= 0.25)
        {
        printf("Duty cycle too high.  Please check.\n");
        abort(1);
        }
      }
    else
      {
     if (((at+mix+trim)/(at+mix+trim+d1)) >= 0.25)
	{
	printf("Duty cycle too high.  Please check.\n");
	abort(1);
	}
      }
   if ((at+mix+trim) >= 0.320)
        { 
        printf(" at+mix+trim > 320msec.  Please check.\n");
        abort(1); 
        }
  }

   if (mult != 0.0)
	{
		if (dmm[C]=='s')
		{
		 printf("dmm[C] must be 'c'.\n");
		 abort(1);
		}
	} 
/* PHASECYCLE CALCULATION */
    ifzero(ssctr);
       hlv(ct,v1); mod2(ct,v3);
    elsenz(ssctr);
       sub(ssval,ssctr,v13);
       hlv(v13,v1); mod2(v13,v3);
    endif(ssctr);

    hlv(v1,v11); hlv(v11,v11);
    add(v1,one,v4);    /* for composite 180 and 1H 180 in direct editting  */
    add(two,v1,v2); 
    dbl(v3,v3); add(v3,v1,v3); assign(v3,oph);
    mod2(v11,v11); dbl(v11,v11);
    assign(v11, v8);	/* for 13C pulse in direct editting */
    add(v11,v1,v11);
    add(one,v11,v12); add(one,v12,v13);
    if (iphase == 2)
	 incr(v3);
    if (iphase == 3)
	 add(v14,v3,v3);             /* TPPI phase increment */

    if ((iphase == 1) || (iphase == 2))
      {
       initval(2.0*(double)((int)(d2*getval("sw1")+0.5)%2),v14);
       add(v3,v14,v3); add(oph,v14,oph);   /*HYPERCOMPLEX with FAD */
      }
	
/* BEGIN ACTUAL PULSE SEQUENCE CODE */
    status(A);
       delay(1.0e-5); rlpower(pwxlvl,DODEV);
       rcvroff();
       if (sstrim != 0.0)
        {
         rlpower((tpwr - 4), TODEV);
	 rgpulse(sstrim,zero,5.0e-5,rof1);
         rgpulse(sstrim,one, 0.0,rof1);
        }
       rlpower(tpwr, TODEV);
       hsdelay(d1);
/* selective saturation period */
       if (satmode[0] == 'y')
         {
          offset(satfrq,TODEV);
          rlpower(satpwr,TODEV);
         rgpulse(satdly,zero,4.0e-5,1.0e-5);
          offset(tof,TODEV);
          rlpower(tpwr,TODEV); rcvroff();
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

          if (satmode[1] == 'y')
            {
             offset(satfrq,TODEV);
             rlpower(satpwr,TODEV);
             rgpulse(null,zero,4.0e-5,1.0e-5); 
             offset(tof,TODEV);
             rlpower(tpwr,TODEV); rcvroff();
             delay(5.0e-5);
            }
          else hsdelay(null);
         }

       add(one,v13,v2);
       rgpulse(pw,v1,0.0,0.0); rcvroff();
       delay(bird); decpulse(pwx,v3);
       delay(d2/2.0);
       rgpulse(2.0*pw,v1,0.0,0.0); rcvroff();
       delay(d2/2.0); decpulse(pwx,v1);
       txphase(v12);
       rlpower(p1lvl,TODEV); 
	if (dm[C] == 'y')
	rlpower(dpwr, DODEV);
       delay(bird - 4.2e-6);
	status(C); 
/* calculate and initialize loop counter */
       if (mix != 0.0)
       {
       cycles = (mix-trim)/(64.66*p1+32*window);
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
            rgpulse(0.66*p1,v12,0.0,0.0);
          endhardloop();
         }
	}
/*	INVERT DIRECT RESPONSES if desired	*/
	if (mult != 0.0)
	{
 	 	rlpower(tpwr,TODEV);
		rlpower(pwxlvl, DODEV);
		delay(bird - 8.4e-6);
		simpulse(2.0*pw,mult*pwx,v4,v8,0.0,0.0); /* orthog. 1H 180 */
		rlpower(dpwr,DODEV);
		delay(bird - 4.2e-6);
	}
	else
	{
		rlpower(dpwr, DODEV);
		delay(1.0e-6);
	}
	status(D);
       txphase(zero); delay(rof2); rcvron();
}
