/*  gcbca_co_nh2.c   with seduce1 decoupling


    This pulse sequence will allow one to perform the following
    experiment:

    3D cbcaconnh correlating cb(i),ca(i) with n(i+1), nh(i+1).
    ref: Grzesiek and Bax, JACS 114 6291, 1992.

(Initially sit at 43 ppm and then jump to 54 ppm - modified by MRB, see below)

    Uses three channels:
         1)  1H             - carrier @ 4.7 ppm [H2O]
         2) 13C             - carrier @ 46 ppm [Cb region] throughout
         3) 15N             - carrier @ 120 ppm  


    Set dm = 'nnn', dmm = 'ccc' 	(modified by MRB)
    Set dm2 = 'nny', dmm2 = 'ccg' 

    Must set phase = 1,2 and phase2 = 1,2 for States-TPPI
    acquisition in t1 [carbon]  and t2 [N].
    
    The flag f1180/f2180 should be set to 'y' if t1 is 
    to be started at halfdwell time. This will give -90, 180
    phasing in f1/f2. If it is set to 'n' the phasing will
    be 0,0 and will still give a perfect baseline.

    Set f1180 to n for (0,0) in C and f2180 = n for (0,0) in N

    Written by Lewis Kay 09/15/92 
    Modified by L.E.K 11/24/92 to include enhanced pfg and shaped C' pulses
    Modified by L.E.K. 01/03/93 to reduce phase cycle by adding gradients
             in the t1 evolution constant time domain. H decoupling now
             begins in the interval immediately following the t1 evolution
             domain.   
    Modified by G.G 10/14/93 to permit DPS and remove gate

    Modified extensively by M.R. Bendall, Varian, March 1995 to include
    better C13 pulses and automated calibration. Uses macro, "c_co_nhcal", 
    and features in:
	   	SLP pulses:     J Magn. Reson. 96, 94-102 (1992)
    		shaka6 composite: Chem. Phys. Lett. 120, 201 (1985) 


                          INSTRUCTIONS FOR USE
         
    1. Calibrate 90 degree pulse lengths for H1 (pw), C13 (pwC) and N15 (pwN)
       at the chosen power levels, tpwr, pwClvl and pwNlvl respectively.

    2. You should check, at least once, the C13 amplifier compression.  At a 
       particular pwClvl (eg 60) and dpwrf=4095, measure the C13 90 - this is
       pwC of course.  Repeat this measurement at dpwrf=2047 (add 6 to dpwr
       if you are using C13 decoupling during this measurement).  A compression
       factor, "comp" is then determined by:
		comp = (C13 90 @ dpwrf=2047) / 2.0 * pwC
       If comp is less than 0.95 serious loss of S/N may occur because the 
       amplitudes of the soft C13 pulses will not be accurately calibrated.
       Either reduce pwClvl until comp > 0.95, or set the value of the comp
       parameter provided to this compression factor. The default is comp=1.0.
       As an example, in one test, setting comp=0.85 as determined above
       reduced the loss in S/N from 30% to 3%. If comp > 0.95, a small
       advantage may result in setting comp to the measured value.

    3. The SLP pulses (shifted laminar pulses - pulses shifted off-resonance)
       must then be calculated using the macro "c_co_nhcal".  These pulses
       are stored in your shapelib and automatically called by the pulse
       sequence code.  This calculation takes into account whether the
       spectrometer is a 500, 600, or 750 system and it only needs to be done
       once for each system.  All other pulse powers and
       pulse lengths are automatically calculated from these settings within
       the pulse sequence code for 500, 600, and 750MHz systems.

    4. The two C13 180 degree pulses after the DIPSI-3 decoupling have been
       replaced by a single 6-pulse composite dubbed "shaka6" which povides
       near perfect refocusing with no phase roll over the Cab region and
       near perfect inversion over the CO region 18kHz off-resonance. 

    5. "c_co_nhcal" also sets ni=ni2=0 and phase=phase2=1 ready for a 1D 
       spectral check.  When ni=0, the calculated lengths of the SLP pulses
       (pwC4, pwC5, pwC7)  are
       displayed by dg to provide confidence to the user that everything has
       been set up properly - do not change these!  Note also that these SLP
       pulse lengths are determined only by dfrq and do not depend on pwC. 

    6. Centre H1 frequency on H2O (4.7ppm), C13 frequency on 46ppm, and N15 
       frequency on the amide region.

    7. The H1 frequency is NOT shifted to the amide region during the sequence
       as suggested in the JMR article as this appears to achieve no purpose.
       The H1 DIPSI2 decoupling was increased to 7.5 kHz to achieve better H2O
       suppression.

    8. The parameter phi7cal (listed in dg when ni=0) is provided to adjust 
       the phase of the last 90 degree C13 pulse, which is phase-shifted by
       the prior 180 degree pulse on the Ca region and by the use of SLP
       pulses on the CO region. The experimentally determined value of this
       phase is also very sensitive to small timing differences (microseconds)
       between the two theta delays.  A value of 267 degrees (or -93) has been
       found for 500 Unity Plus,  but this will change significantly
       with other systems.  Check this phase via 1D spectra - maximise signal,
       or for greater accuracy obtain the value for a null and then add or 
       subtract 90 degrees.

    9. Most gradient levels and durations require no adjustment.  The coherence
       pathway gradients gzlvl1 and gzlvl2 should be adjusted for maximum
       signal in a 1D spectrum (normally gzlvl2 is slightly less than gzlvl1).
       typical values are 15000/14850 with the duration gt1=0.0025.

   10. All echo times within the sequence have been adjusted to provide maximum
       signal for the protein alphalytic protease (mw 2?000).  These are not
       under parameter adjustment, but can be changed within the sequence code.
       To do this replace each pair in turn with d4 and maximise 1D signals with
       an array of d4.
*/

#include <standard.h>
#include <math.h>
#define PI 3.1416

static int  phi1[1]  = {0},
            phi2[1]  = {1},
            phi3[1]  = {0},
            phi4[1]  = {0},
            phi5[2]  = {0,2},
            phi6[2]  = {2,0},
            phi7[1]  = {0},
            phi8[1] = {0},
            phi9[4] = {0,0,2,2},
            rec[2] = {0,2},
            phi11[1] = {0},
            phi12[1] = {1},
            phi13[1] = {1};
           
static double d2_init=0.0, d3_init=0.0;
            
pulsesequence()
{
/* DECLARE VARIABLES */

 char    /*   fsat[MAXSTR],
	    fscuba[MAXSTR],  */
            f1180[MAXSTR],    /* Flag to start t1 @ halfdwell             */
            f2180[MAXSTR];    /* Flag to start t2 @ halfdwell             */
 
 int         phase, phase2, ni, ni2, icosel,  /* used to get n and p type */
             t1_counter,   /* used for states tppi in t1           */ 
             t2_counter;   /* used for states tppi in t2           */ 

 double      tau1,         /*  t1 delay */
             tau2,         /*  t2 delay */
  /*         tsatpwr,       low level 1H trans.power for presat  */
            
 	    pwClvl,	     /* coarse power for C13 pulses		    */
            pwC,             /* C13 90 degree pulse length at pwClvl        */
	    rf0,	     /* maximum fine power when using pwC pulses    */

/* 90 degree pulse at Cab(46ppm), first off-resonance null at CO (177ppm)   */
            pwC1,            /* 90 degree pulse length on C13 at rf1        */
            rf1,             /* fine power for 5.1 kHz rf for 600MHz magnet */

/* 180 degree pulse at Cab(46ppm), first off-resonance null at CO(177ppm)   */
            pwC2,            /* 180 degree pulse length at rf2              */
            rf2,             /* fine power for 11.4 kHz rf for 600MHz magnet*/

/* p_d is used to calculate the selective decoupling on the Cab region      */
            p_d,             /* 50 degree pulse for DIPSI-3 at rf3 or dpwr1 */
            rf3,             /* fine power for 6.75 kHz rf for 600MHz magnet*/
	    pwCd,	     /* 90 degree pulse length of seduce1 element   */
	    sedpwr,          /* fine power for seduce1 decoupling on Ca   */

/* the following pulse lengths for SLP pulses are automatically calculated  */
/* by the macro "c_co_nhcal".  SLP pulse shapes, "offC4" etc are called     */
/* directly from your shapelib.                    			    */
            pwC4,       /* 180 degree pulse at Ca(56ppm) null at CO(177ppm) */
	    pwC5,	/* 90 degree pulse at CO(177ppm) null at Ca(56 ppm) */
	    pwC7,	/* 180 degree selective sinc pulse on CO(177ppm)    */
	    rf4,	/* fine power for the pwC4 ("offC4") pulse	    */
	    rf5,	/* fine power for the pwC5 ("offC5") pulse	    */	
	    rf7,	/* fine power for the pwC7 ("offC7") pulse	    */

	    comp,      	     /* adjustment for C13 amplifier compression    */

	    pwH,	     /* H1 90 degree pulse length at tpwr1	  */
	    tpwr1,	     /* 5 kHz rf for 600MHz magnet for DIPSI-2    */
	    DIPSI2time,      /* total length of DIPSI-2 decoupling	  */

	    pwNlvl,	     /* power for N15 pulses			  */
            pwN,             /* N15 90 degree pulse length at pwNlvl      */

	    phi7cal,        /* phase in degrees of the last C13 90 pulse  */

	    sw1,
	    sw2,

	    gt1, gzlvl1, gzlvl2;  /* coherence pathway gradients          */

/* LOAD VARIABLES */
  getstr("f1180",f1180);
  getstr("f2180",f2180);

/*  getstr("fsat",fsat);*/
/*  getstr("fscuba",fscuba);*/
/*  tsatpwr = getval("tsatpwr");*/

	pwClvl = getval("pwClvl");
	pwC = getval("pwC");
	pwC4 = getval("pwC4");
	pwC5 = getval("pwC5");
	pwC7 = getval("pwC7");
	comp = getval("comp");
	pwNlvl = getval("pwNlvl");
    	pwN = getval("pwN");
        ni = getval("ni");
	ni2 = getval("ni2");
	sw1 = getval("sw1");
	sw2 = getval("sw2");
	phi7cal = getval("phi7cal");
	gt1 = getval("gt1");
	gzlvl1 = getval("gzlvl1");
	gzlvl2 = getval("gzlvl2");

	phase = (int) (getval("phase") + 0.5);
	phase2 = (int) (getval("phase2") + 0.5);

/* LOAD PHASE TABLE */

  settable(t1,1,phi1);
  settable(t2,1,phi2);
  settable(t3,1,phi3);
  settable(t4,1,phi4);
  settable(t5,2,phi5);
  settable(t6,2,phi6);
  settable(t7,1,phi7);
  settable(t8,1,phi8);
  settable(t9,4,phi9);
  settable(t10,2,rec);
  settable(t11,1,phi11);
  settable(t12,1,phi12);
  settable(t13,1,phi13);

/*   INITIALIZE VARIABLES   */

       	if( dpwrf < 4095 )
	{
		printf("reset dpwrf=4095 and recalibrate C13 90 degree pulse");
		abort(1);
	}

	/* maximum fine power for pwC pulses */
	rf0 = 4095.0;

	/* 90 degree pulse on Cab, null at CO 131ppm away */
	pwC1 = sqrt(15.0)/(4.0*131.0*dfrq);
        rf1 = (comp*4095.0*pwC)/pwC1;
	rf1 = (int) (rf1 + 0.5);
	
	/* 180 degree pulse on Cab, null at CO 131ppm away */
        pwC2 = sqrt(3.0)/(2.0*131.0*dfrq);
	rf2 = (4095.0*pwC*2.0)/pwC2;
	rf2 = (int) (rf2 + 0.5);	
	if( rf2 > 4095 )
	{
		printf("retune so that C13 90 is less than 21.9us*600/sfrq");
		abort(1);
	}

	/* dipsi-3 decoupling */	
 	p_d = (20.6e-6*600.0)/sfrq;
 	rf3 = (comp*4095.0*pwC*5.0)/(p_d*9.0);
	rf3 = (int) (rf3 + 0.5);
	
	/* 180 degree pulse on Ca, null at CO 131ppm away */
	rf4 = (comp*4095.0*pwC*2.0)/pwC4;
	rf4 = (int) (rf4 + 0.5);

	/* 90 degree one-lobe sinc pulse on CO, null at Cab 131ppm away */
	rf5 = (comp*4095.0*pwC*1.65)/pwC5;	/* needs 1.65 times more     */
	rf5 = (int) (rf5 + 0.5);		/* power than a square pulse */

	/* 180 degree one-lobe sinc pulse on CO, null at Cab 131ppm away */
	rf7 = (comp*4095.0*pwC*2.0*1.65)/pwC7;	/* needs 1.65 times more     */
	rf7 = (int) (rf7 + 0.5);		/* power than a square pulse */

	pwCd = (275.0e-6*600.0)/sfrq;
	sedpwr =(comp*4095.0*pwC)/(pwCd*0.47);
	sedpwr = (int) (sedpwr + 0.5); 	

	DIPSI2time = 2.0*3.0e-3 + 2.0*14.4e-3 + 2.0*14.0e-3 - 5.4e-3 + 0.5*pwC1 + pwC5 + pwN + 2*1.0638e-3;
	pwH = (DIPSI2time*90.0)/(12.0*2590.0*4.0);   /* 12 cycles of DIPSI-2  */
	tpwr1 = tpwr - 20.0*log10(pwH/pw);
	tpwr1 = (int) (tpwr1 + 0.5);
 
/* CHECK VALIDITY OF PARAMETER RANGES */

    if( 0.5*ni*1/(sw1) > 3.3e-3 - 0.25e-3 - 24.15e-6 - pwC7 )
    {
        printf(" ni is too big\n");
        abort(1);
    }

    if( ni2*1/(sw2) > 2.0*11.0e-3)
    {
        printf(" ni2 is too big\n");
        abort(1);
    }

    if((dm[A] == 'y' || dm[B] == 'y' || dm[C] == 'y' ))
    {
        printf("incorrect dec1 decoupler flags! Should be 'nnn' ");
        abort(1);
    }

    if((dm2[A] == 'y' || dm2[B] == 'y'))
    {
        printf("incorrect dec2 decoupler flags! Should be 'nny' ");
        abort(1);
    }


/*    if( tsatpwr > 6 )
    {
        printf("TSATPWR too large !!!  ");
        abort(1);
    }*/


    if( dpwr2 > 46 )
    {
        printf("don't fry the probe, DPWR2 too large!  ");
        abort(1);
    }

    if( pw > 20.0e-6 )
    {
        printf("dont fry the probe, pw too high ! ");
        abort(1);
    } 
  
    if( pwN > 100.0e-6 )
    {
        printf("dont fry the probe, pwN too high ! ");
        abort(1);
    } 


/*  Phase incrementation for hypercomplex 2D data */

    if (phase == 2)
      tsadd(t3,1,4);  
    if (phase2 == 2) {
      tsadd(t11,2,4);   
      icosel = 1;
    }
    else icosel = -1; 

/*  Set up f1180  tau1 = t1               */
   
    tau1 = d2;
    if(f1180[A] == 'y') {
        tau1 += ( 1.0 / (2.0*sw1) );
        if(tau1 < 0.2e-6) tau1 = 0.0;
    }
        tau1 = tau1/2.0;

/*  Set up f2180  tau2 = t2               */

    tau2 = d3;
    if(f2180[A] == 'y') {
        tau2 += ( 1.0 / (2.0*sw2) ); 
        if(tau2 < 0.2e-6) tau2 = 0.0;
    }
        tau2 = tau2/2.0;

/* Calculate modifications to phases for States-TPPI acquisition          */

   if( ix == 1) d2_init = d2 ;
   t1_counter = (int) ( (d2-d2_init)*sw1 + 0.5 );
   if(t1_counter % 2) {
      tsadd(t3,2,4);     
      tsadd(t10,2,4);    
    }

   if( ix == 1) d3_init = d3 ;
   t2_counter = (int) ( (d3-d3_init)*sw2 + 0.5 );
   if(t2_counter % 2) {
      tsadd(t8,2,4);  
      tsadd(t10,2,4);    
    }

/* BEGIN ACTUAL PULSE SEQUENCE */

status(A);
  rcvroff();


/* Presaturation Period - commented out 
   if (fsat[A] == 'y')
   {
        rlpower(tsatpwr,TODEV);  Set transmitter power for 1H presaturation 
	delay(2.0e-5);
        rgpulse(d1,zero,rof1,rof1);
   	rlpower(tpwr,TODEV);       Set transmitter power for hard 1H pulses 
	delay(2.0e-5);
	if(fscuba[A] == 'y')
	{
		delay(2.2e-2);
		rgpulse(pw,zero,2.0e-6,0.0);
		rgpulse(2*pw,one,2.0e-6,0.0);
		rgpulse(pw,zero,2.0e-6,0.0);
		delay(2.2e-2);
	}
   }
   else*/
          
   	delay(d1);

	rlpower(tpwr, TODEV);
	rlpower(pwClvl, DODEV);
 	rlpower(pwNlvl, DO2DEV);
	rlpwrf(rf0, DODEV);
	offset(tof, TODEV);
	txphase(zero);
   	delay(1.0e-5);

/* Begin Pulses */

status(B);

/* ensure that magnetization originates on H and not 13C  */
   	decrgpulse(pwC,zero,0.0,0.0);
   	delay(0.2e-6);
   	rgradient('z',8000.0);
   	delay(0.5e-3);
   	rgradient('z',0.0);
   	delay(150.0e-6);
   	decrgpulse(pwC,one,0.0,0.0);
   	delay(0.2e-6);
   	rgradient('z',8000.0);
   	delay(0.5e-3);
   	rgradient('z',0.0);
   	delay(150.0e-6);

   	rgpulse(pw,zero,0.0,0.0);                    /* 90 deg 1H pulse */
   	decphase(t1);

   	rgradient('z',8000.0);					/*  6.9us  */
   	delay(0.5e-3);
   	rgradient('z',0.0);					/*  6.9us  */

   	delay(1.7e-3 - 0.5*pw - 0.5138e-3);
                       
   	simpulse(2.0*pw,2.0*pwC,zero,t1,0.0,0.0);

   	rgradient('z',8000.0);					/*  6.9us  */
   	delay(0.5e-3);
   	rgradient('z',0.0);					/*  6.9us  */

   	txphase(t2);
 	decphase(t3);
 
   	delay(1.7e-3 - 0.5*pwC - 0.5138e-3);

   	rgpulse(pw,t2,0.0,0.0);
   	rgradient('z',20000.0);
   	delay(1.0e-3);
   	rgradient('z',0.0);
   	delay(50.0e-6);
   	decrgpulse(pwC,t3,0.0,0.0);				/*  point a  */

     	delay(tau1);

        rlpwrf(rf7, DODEV);					/*  4.6us    */
								/*  17.25us  */
	sim3shaped_pulse("square", "offC7", "square", 0.0, pwC7, 2.0*pwN, zero, zero, zero, 0.0, 0.0);

     	rgradient('z',20000.0);					/*  6.9us  */
     	delay(0.25e-3);
     	rgradient('z',0.0);					/*  6.9us  */

        rlpwrf(rf2, DODEV);					/*  4.6us    */

     	delay(1.1e-3 - 2.0*pwN - 0.25e-3 - 40.25e-6);
     
     	rgpulse(2.0*pw,zero,0.0,0.0);

     	decphase(t4);
     	delay(3.3e-3 - 1.1e-3);

     	decrgpulse(pwC2,t4,0.0,0.0);

     	rgradient('z',20000.0);					/*  6.9us  */
     	delay(0.25e-3);
     	rgradient('z',0.0);					/*  6.9us  */

     	delay(3.3e-3 - tau1 - 0.25e-3 - 28.75e-6 - pwC7);

        rlpwrf(rf7, DODEV);					/*  4.6us    */
 								/*  5.75us   */
     	decshaped_pulse("offC7",pwC7,zero,0.0,0.0);
								/*  point b  */
	rlpwrf(rf1, DODEV);					/*  4.6us    */
   	decrgpulse(pwC1,zero,2.0e-6,0.0);
								/*  point c  */

	rlpower(tpwr1, TODEV);					/*  2.3us    */
	obsprgon("dipsi2", pwH, 5.0);				/*  5.75us   */
	xmtron();

	rlpwrf(rf0, DODEV);
	decphase(t5);
	delay(3.0e-3 - 12.65e-6 - 0.5*10.933*pwC);

	decrgpulse(pwC*158.0/90.0, t5, 0.0, 0.0);
	decrgpulse(pwC*171.2/90.0, t6, 0.0, 0.0);
	decrgpulse(pwC*342.8/90.0, t5, 0.0, 0.0);	/* Shaka composite   */
	decrgpulse(pwC*145.5/90.0, t6, 0.0, 0.0);
	decrgpulse(pwC*81.2/90.0, t5, 0.0, 0.0);
	decrgpulse(pwC*85.3/90.0, t6, 0.0, 0.0);

	rlpwrf(rf1, DODEV);					/*  4.6us    */
	decphase(zero);
	delay(3.0e-3 - 4.6e-6 - 0.5*10.933*pwC - 0.5*pwC1);
								/*  point d  */
	decrgpulse(pwC1, zero, 0.0, 0.0);
	decphase(t5);
	rgradient('z', 20000.0);				/*  6.9us    */
	delay(1.0e-3);
	rgradient('z', 0.0);					/*  6.9us    */
	delay(5.0e-5);
        rlpwrf(rf5, DODEV);					/*  4.6us    */
								/*  5.75us   */
	decshaped_pulse("offC5", pwC5, t5, 0.0, 0.0);	
 								/*  point e  */
        rlpwrf(rf4, DODEV);					/*  4.6us    */
 	decphase(zero);
	delay(4.6e-3);
								/*  5.75us   */
	decshaped_pulse("offC4", pwC4, zero, 0.0, 0.0);

	dec2phase(zero);
        rlpwrf(rf7, DODEV);					/*  4.6us    */

	delay(14.4e-3 - 4.6e-3 - 32.2e-6 - 0.5*pwC5 - pwC4 - pwN);
								/*  17.25us  */
	sim3shaped_pulse("square", "offC7", "square", 0.0, pwC7, 2.0*pwN, zero, zero, zero, 0.0, 0.0); 

        rlpwrf(rf5, DODEV);					/*  4.6us    */
	initval(phi7cal, v7);
	stepsize(1.0, DODEV);
	dcplrphase(v7);						/*  1.15us   */
	dec2phase(t8);

	delay(14.4e-3 - 11.5e-6 - pwN - 0.5*pwC5);
								/*  point f  */
								/*  5.75us   */
	decshaped_pulse("offC5", pwC5, t7, 0.0, 0.0);
	rgradient('z', 20000.0);				/*  6.9us    */
	delay(1.0e-3);
	rgradient('z', 0.0);					/*  6.9us    */
	delay(5.0e-5);
        rlpwrf(rf7, DODEV);					/*  4.6us    */
	dec2rgpulse(pwN, t8, 0.0, 0.0);

	dcplrphase(zero);					/*  1.15us   */
	decphase(zero);
	dec2phase(t9);

	rlpwrf(sedpwr, DODEV);					/*  4.6us    */
	decprgon("seduce1", pwCd, 0.9);				/*  5.75us   */
	decon();

	delay(14.0e-3 - 33.35e-6 - 1.5*pwN - tau2);

	decoff();
	decprgoff();
	rlpwrf(rf7, DODEV);					/*  4.6us    */
								/* 17.25     */	sim3shaped_pulse("square", "offC7", "square", 0.0, pwC7, 2.0*pwN, zero, zero, t9, 0.0, 0.0);

	dec2phase(t11);
	txphase(zero);

	rlpwrf(sedpwr, DODEV);					/*  4.6us    */
	decprgon("seduce1", pwCd, 0.9);				/*  5.75us   */
	decon();

	delay(14.0e-3 - pwN - 10.35e-6 - 5.4e-3 + tau2);

	xmtroff();
	obsprgoff();

	delay(5.4e-3 - 16.1e-6 - gt1 - 1.0e-4 - 0.5*pwN);
	decoff();
	decprgoff();

	rgradient('z', icosel*gzlvl1);				/*  6.9us    */
	delay(gt1);
	rgradient('z', 0.0);					/*  6.9us    */
	delay(1.0e-4);

	rlpower(tpwr, TODEV);					/*  2.3us    */
								/*  point g  */
	sim3pulse(pw, 0.0, pwN, zero, zero, t11, 0.0, 0.0);

	dec2phase(zero);

	rgradient('z', 2000.0);					/*  6.9us    */
	delay(0.5e-3);
	rgradient('z', 0.0);					/*  6.9us    */

	delay(2.4e-3 - 1.5*pwN - 0.5138e-3);

	sim3pulse(2.0*pw, 0.0, 2.0*pwN, zero, zero, zero, 0.0, 0.0);

	rgradient('z', 2000.0);					/*  6.9us    */
	delay(0.5e-3);
	rgradient('z', 0.0);					/*  6.9us    */

	txphase(t12);
	dec2phase(t13);
	delay(2.4e-3  - 1.5*pwN - 0.5138e-3);

	sim3pulse(pw, 0.0, pwN, t12, zero, t13, 0.0, 0.0);

	txphase(zero);
	dec2phase(zero);

	rgradient('z', 3000.0);					/*  6.9us    */
	delay(0.5e-3);
	rgradient('z', 0.0);					/*  6.9us    */

	delay(2.4e-3 - 1.5*pwN - 0.5138e-3);
	
	sim3pulse(2.0*pw, 0.0, 2.0*pwN, zero, zero, zero, 0.0, 0.0);

	rgradient('z', 3000.0);					/*  6.9us    */
	delay(0.5e-3);
	rgradient('z', 0.0);					/*  6.9us    */

	delay(2.4e-3 - pwN - 0.5*pw - 0.5138e-3);

	rgpulse(pw, zero, 0.0, 0.0);

	delay((gt1/10.0) + 5.0e-5 - 0.5*pw + 16.1e-6);
	
	rgpulse(2.0*pw, zero, 0.0, 0.0);

	rlpower(dpwr2, DO2DEV);					/*  2.3us    */

				      	
	rgradient('z', gzlvl2);					/*  6.9us    */
	delay(gt1/10.0);
	rgradient('z', 0.0);					/*  6.9us    */

	delay(5.0e-5);		

status(C);
	setreceiver(t10);

}		 

