/*  c_co_nh.c                

    pulse sequence: J Magn. Reson.  B  101, 114-119 (1993)
    SLP pulses:     J Magn. Reson. 96, 94-102 (1992)
    shaka6 composite:  Chem. Phys. Lett. 120, 201 (1985) 

    written by Robin Bendall, Varian, March 95  (version 20 something)


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
       frequency on the amide region.  The carrier remains at 46ppm, ie at 
       Cb throughout the sequence.

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

    9. All echo times within the sequence have been adjusted to provide maximum
       signal for the protein alphalytic protease (mw 2?000).  These are not
       under parameter adjustment, but can be changed within the sequence code.
       To do this replace each pair in turn with d4 and maximise 1D signals with
       an array of d4. 

  
  
*/

#include <standard.h>
#include <math.h>



static int   phi1[1]  = {1},
             phi2[2]  = {0,2},
             phi3[1]  = {0},
             phi4[32] = {0,0,0,0,0,0,0,0,1,1,1,1,1,1,1,1,
                         2,2,2,2,2,2,2,2,3,3,3,3,3,3,3,3},
             phi5[8]  = {0,0,0,0,2,2,2,2},
             phi6[4]  = {0,0,2,2},  
        /* phi7 of 48.5 (recalibrated to 267 degrees) degrees set by phi7cal */
             phi8[8]  = {2,2,2,2,0,0,0,0},
             phi9[16] = {0,0,0,0,0,0,0,0,2,2,2,2,2,2,2,2},
             phi10[32]= {0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
                         2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2}, 
             rec[32]  = {0,2,2,0,2,0,0,2,2,0,0,2,0,2,2,0,
                         2,0,0,2,0,2,2,0,0,2,2,0,2,0,0,2};

static double   d2_init=0.0, d3_init=0.0;



pulsesequence()
{


/*  DECLARE VARIABLES  */

char     fad1[MAXSTR], fad2[MAXSTR];

int   	 phase, phase2,
      	 t1_counter, t2_counter,
	 ni;

           
           
double      pwClvl,	     /* coarse power for C13 pulses		    */
            pwC,             /* C13 90 degree pulse length at pwClvl        */
	    rf0,	     /* maximum fine power when using pwC pulses    */

/* 90 degree pulse at Cab(46ppm), first off-resonance null at CO (177ppm)   */
            pwC1,            /* 90 degree pulse length on C13 at rf1        */
            rf1,             /* fine power for 5.1 kHz rf for 600MHz magnet */

/* 180 degree pulse at Cab(46ppm), first off-resonance null at CO(177ppm)   */
            pwC2,            /* 180 degree pulse length at rf2              */
            rf2,             /* fine power for 11.4 kHz rf for 600MHz magnet*/

/* p_d is used to calculate the selective decoupling on the Cab region      */
            p_d,             /* 50 degree pulse for DIPSI-3 at rf3          */
            rf3,             /* fine power for 6.75 kHz rf for 600MHz magnet*/

/* the following pulse lengths for SLP pulses are automatically calculated  */
/* by the macro "c_co_nhcal".  SLP pulse shapes, "offC4" etc are called     */
/* directly from your shapelib.                    			    */
            pwC4,       /* 180 degree pulse at Ca(56ppm) null at CO(177ppm) */
	    pwC5,	/* 90 degree pulse at CO(177ppm) null at Ca(56 ppm) */
	    pwC7,	/* 180 degree selective sinc pulse on CO(177ppm)    */
	    rf4,	/* fine power for the pwC4 ("offC4") pulse	    */
	    rf5,	/* fine power for the pwC5 ("offC5") pulse	    */	
	    rf7,	/* fine power for the pwC7 ("offC7") pulse	    */

	    comp,	/* adjustment for C13 amplifier compression	    */

	    pwH,	     /* H1 90 degree pulse length at tpwr1	  */
	    tpwr1,	     /* 5 kHz rf for 600MHz magnet for DIPSI-2    */
	    DIPSI2time,      /* total length of DIPSI-2 decoupling	  */

	    pwNlvl,	     /* power for N15 pulses			  */
            pwN,             /* N15 90 degree pulse length at pwNlvl      */

            t1a,             /* time increments for first dimension       */
	    t1b,
	    t1c,

	    phi7cal,        /* phase in degrees of the last C13 90 pulse  */

	    sw1,
	    sw2;
 
            

/*  LOAD VARIABLES  */

	pwClvl = getval("pwClvl");
	pwC = getval("pwC");
	pwC4 = getval("pwC4");
	pwC5 = getval("pwC5");
	pwC7 = getval("pwC7");
	comp = getval("comp");
	pwNlvl = getval("pwNlvl");
    	pwN = getval("pwN");
        ni = getval("ni");
	sw1 = getval("sw1");
	sw2 = getval("sw2");
	phi7cal = getval("phi7cal");

	phase = (int) (getval("phase") + 0.5);
	phase2 = (int) (getval("phase2") + 0.5);

	getstr("fad1", fad1);
	getstr("fad2", fad2);

        

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

	DIPSI2time = 2.0*3.0e-3 + 2.0*14.4e-3 + 2.0*14.0e-3 - 5.4e-3 + 0.5*pwC1 + 0.5*pwC5 + 0.5*pwN;
	pwH = (DIPSI2time*90.0)/(12.0*2590.0*4.0);   /* 12 cycles of DIPSI-2  */
	tpwr1 = tpwr - 20.0*log10(pwH/pw);
	tpwr1 = (int) (tpwr1 + 0.5);





/*   LOAD PHASE TABLE    */


	settable(t1,1,phi1);
	settable(t2,2,phi2);
	settable(t3,1,phi3);
	settable(t4,32,phi4);
	settable(t5,8,phi5);
	settable(t6,4,phi6);

	settable(t8,8,phi8);
	settable(t9,16,phi9);
	settable(t10,32,phi10);
	settable(t11,32,rec);


/*   C13 TIME INCREMENTATION CALCULATIONS AND PHASE INCREMENTATIONS  */

	if( ix == 1)   d2_init = d2;
        t1_counter = (int) ((d2 - d2_init)*sw1 + 0.5);
         
	t1a = d2/2.0;
	t1b = (d2/2.0) - 0.95e-3*t1_counter/(ni - 1);
     	t1c = 0.95e-3 -  0.95e-3*t1_counter/(ni - 1);


	/* Add in States-Haberkorn element */

	if (phase == 2)    tsadd(t3,1,4);
	if (phase2 == 2)   tsadd(t8,1,4);	  


	/* Add in FAD */

	if (fad1[A] == 'y')
	{
	if (ix == 1)  d2_init = d2;
	t1_counter = (int) ((d2 - d2_init)*sw1 + 0.5);
	if (t1_counter % 2)
		{
		tsadd(t3,2,4);
		tsadd(t11,2,4);
		}
	}

	if (fad2[A] == 'y')
	{
	if (ix == 1)   d3_init = d3;
	t2_counter = (int) ((d3 - d3_init)*sw2 + 0.5);
	if (t2_counter % 2)
		{
		tsadd(t8,2,4);
		tsadd(t11,2,4);
		}
	}



/*   BEGIN PULSE SEQUENCE   */

status(A);
	hsdelay(d1);

status(B);
	rcvroff();
	rlpower(tpwr, TODEV);
	rlpower(pwClvl, DODEV);
 	rlpower(pwNlvl, DO2DEV);
	rlpwrf(rf0, DODEV);
	offset(tof, TODEV);
	txphase(t1);
	delay(2.5e-5);

	rgpulse(pw, t1, 0.0, 0.0);
								/* point a */
        txphase(zero);
        decphase(t2);
	delay(1.7e-3 - 0.5*pw);

        simpulse(2*pw, 2.0*pwC, zero, t2, 0.0, 0.0);

        txphase(t1);
	decphase(t3);
	delay(1.7e-3);

	rgpulse(1.8e-3, t1, 0.0, 0.0);				/* point b */
        simpulse(pw, pwC, t2, t3, 0.3e-6, 0.0);
								/* point c */
	decphase(zero);
	delay(t1a);

	rlpwrf(rf7, DODEV);
	decshaped_pulse("offC7", pwC7, zero, 0.0, 0.0);
        rlpwrf(rf2, DODEV);

	txphase(zero);
	delay(0.95e-3 - 2.0*pw);

	rgpulse(2.0*pw, zero, 0.0, 0.0);

	decphase(t4);
	delay(t1b);

	decrgpulse(pwC2, t4, 0.0, 0.0);

	decphase(zero);
	delay(t1c);

	rlpwrf(rf7, DODEV);
	decshaped_pulse("offC7", pwC7, zero, 0.0, 0.0);
	rlpwrf(rf3, DODEV);
								/* point d */
	decrgpulse(1.0e-3, zero, 0.0, 0.0);

	initval(2.0, v2);
	starthardloop(v2);

     decrgpulse(4.9*p_d,one,0.0,0.0);
     decrgpulse(7.9*p_d,three,0.0,0.0);
     decrgpulse(5.0*p_d,one,0.0,0.0);
     decrgpulse(5.5*p_d,three,0.0,0.0);
     decrgpulse(0.6*p_d,one,0.0,0.0);
     decrgpulse(4.6*p_d,three,0.0,0.0);
     decrgpulse(7.2*p_d,one,0.0,0.0);
     decrgpulse(4.9*p_d,three,0.0,0.0);
     decrgpulse(7.4*p_d,one,0.0,0.0);
     decrgpulse(6.8*p_d,three,0.0,0.0);
     decrgpulse(7.0*p_d,one,0.0,0.0);
     decrgpulse(5.2*p_d,three,0.0,0.0);
     decrgpulse(5.4*p_d,one,0.0,0.0);
     decrgpulse(0.6*p_d,three,0.0,0.0);
     decrgpulse(4.5*p_d,one,0.0,0.0);
     decrgpulse(7.3*p_d,three,0.0,0.0);
     decrgpulse(5.1*p_d,one,0.0,0.0);
     decrgpulse(7.9*p_d,three,0.0,0.0);

     decrgpulse(4.9*p_d,three,0.0,0.0);
     decrgpulse(7.9*p_d,one,0.0,0.0);
     decrgpulse(5.0*p_d,three,0.0,0.0);
     decrgpulse(5.5*p_d,one,0.0,0.0);
     decrgpulse(0.6*p_d,three,0.0,0.0);
     decrgpulse(4.6*p_d,one,0.0,0.0);
     decrgpulse(7.2*p_d,three,0.0,0.0);
     decrgpulse(4.9*p_d,one,0.0,0.0);
     decrgpulse(7.4*p_d,three,0.0,0.0);
     decrgpulse(6.8*p_d,one,0.0,0.0);
     decrgpulse(7.0*p_d,three,0.0,0.0);
     decrgpulse(5.2*p_d,one,0.0,0.0);
     decrgpulse(5.4*p_d,three,0.0,0.0);
     decrgpulse(0.6*p_d,one,0.0,0.0);
     decrgpulse(4.5*p_d,three,0.0,0.0);
     decrgpulse(7.3*p_d,one,0.0,0.0);
     decrgpulse(5.1*p_d,three,0.0,0.0);
     decrgpulse(7.9*p_d,one,0.0,0.0);

     decrgpulse(4.9*p_d,three,0.0,0.0);
     decrgpulse(7.9*p_d,one,0.0,0.0);
     decrgpulse(5.0*p_d,three,0.0,0.0);
     decrgpulse(5.5*p_d,one,0.0,0.0);
     decrgpulse(0.6*p_d,three,0.0,0.0);
     decrgpulse(4.6*p_d,one,0.0,0.0);
     decrgpulse(7.2*p_d,three,0.0,0.0);
     decrgpulse(4.9*p_d,one,0.0,0.0);
     decrgpulse(7.4*p_d,three,0.0,0.0);
     decrgpulse(6.8*p_d,one,0.0,0.0);
     decrgpulse(7.0*p_d,three,0.0,0.0);
     decrgpulse(5.2*p_d,one,0.0,0.0);
     decrgpulse(5.4*p_d,three,0.0,0.0);
     decrgpulse(0.6*p_d,one,0.0,0.0);
     decrgpulse(4.5*p_d,three,0.0,0.0);
     decrgpulse(7.3*p_d,one,0.0,0.0);
     decrgpulse(5.1*p_d,three,0.0,0.0);
     decrgpulse(7.9*p_d,one,0.0,0.0);

     decrgpulse(4.9*p_d,one,0.0,0.0);
     decrgpulse(7.9*p_d,three,0.0,0.0);
     decrgpulse(5.0*p_d,one,0.0,0.0);
     decrgpulse(5.5*p_d,three,0.0,0.0);
     decrgpulse(0.6*p_d,one,0.0,0.0);
     decrgpulse(4.6*p_d,three,0.0,0.0);
     decrgpulse(7.2*p_d,one,0.0,0.0);
     decrgpulse(4.9*p_d,three,0.0,0.0);
     decrgpulse(7.4*p_d,one,0.0,0.0);
     decrgpulse(6.8*p_d,three,0.0,0.0);
     decrgpulse(7.0*p_d,one,0.0,0.0);
     decrgpulse(5.2*p_d,three,0.0,0.0);
     decrgpulse(5.4*p_d,one,0.0,0.0);
     decrgpulse(0.6*p_d,three,0.0,0.0);
     decrgpulse(4.5*p_d,one,0.0,0.0);
     decrgpulse(7.3*p_d,three,0.0,0.0);
     decrgpulse(5.1*p_d,one,0.0,0.0);
     decrgpulse(7.9*p_d,three,0.0,0.0);

	endhardloop();						/*  point e  */

	rlpower(tpwr1, TODEV);					/*  2.3us    */
	obsprgon("dipsi2", pwH, 5.0);				/*  5.75us   */
	xmtron();

	rlpwrf(rf0, DODEV);					/*  4.6us    */
	decphase(t5);
	delay(3.0e-3 - 12.65e-6 - 0.5*10.933*pwC);

	decrgpulse(pwC*158.0/90.0, t5, 0.0, 0.0);
	decrgpulse(pwC*171.2/90.0, t8, 0.0, 0.0);
	decrgpulse(pwC*342.8/90.0, t5, 0.0, 0.0);	/* Shaka composite   */
	decrgpulse(pwC*145.5/90.0, t8, 0.0, 0.0);
	decrgpulse(pwC*81.2/90.0, t5, 0.0, 0.0);
	decrgpulse(pwC*85.3/90.0, t8, 0.0, 0.0);

	rlpwrf(rf1, DODEV);					/*  4.6us    */
	decphase(zero);
	delay(3.0e-3 - 4.6e-6 - 0.5*10.933*pwC - 0.5*pwC1);
								/*  point f  */
	decrgpulse(pwC1, zero, 0.0, 0.0);
	decphase(t6);
	rlpwrf(rf5, DODEV);					/*  4.6us    */
								/*  5.75us   */
	decshaped_pulse("offC5", pwC5, t6, 0.0, 0.0);	
 								/*  point g  */
	rlpwrf(rf4, DODEV);					/*  4.6us    */
	decphase(zero);
	delay(4.6e-3);
								/*  5.75us   */
	decshaped_pulse("offC4", pwC4, zero, 0.0, 0.0);

	rlpwrf(rf7, DODEV);					/*  4.6us    */
	dec2phase(zero);
	delay(14.4e-3 - 4.6e-3 - 32.2e-6 - 0.5*pwC5 - pwC4 - pwN);
								/*  17.25us  */
	sim3shaped_pulse("square", "offC7", "square", 0.0, pwC7, 2.0*pwN, zero, zero, zero, 0.0, 0.0);
 
	rlpwrf(rf5, DODEV);					/*  4.6us    */
	initval(phi7cal, v7);
	stepsize(1.0, DODEV);
	dcplrphase(v7);						/*  1.15us   */
	dec2phase(t5);

	delay(14.4e-3 - 23.0e-6 - pwN - 0.5*pwC5);
								/*  point h  */
								/*  17.25us  */
	sim3shaped_pulse("square", "offC5", "square", 0.0, pwC5, pwN, zero, zero, t5, 0.0, 0.0);
	dcplrphase(zero);					/*  1.15us   */
								/*  point i  */
	rlpwrf(rf7, DODEV);					/*  4.6us    */
	decphase(zero);
	dec2phase(t9);

	delay(14.0e-3 - 20.5e-6 - 0.5*pwC5 + 0.5*pwN - pwN - (d3/2));
								/* 17.25us  */
	sim3shaped_pulse("square", "offC7", "square", 0.0, pwC7, 2.0*pwN, zero, zero, t9, 0.0, 0.0);

	rlpwrf(rf4, DODEV);					/*  4.6us    */
	dec2phase(t10);
	txphase(zero);

    if ( (d3/2) > 5.4e-3 )
		{
		delay(14.0e-3 - pwN - 4.6e-6 - pwC4 - 5.75e-6);
								/*  5.75us   */
		decshaped_pulse("offC4", pwC4, zero, 0.0, 0.0);
		delay((d3/2) - 5.4e-3);
		xmtroff();
		obsprgoff();
		delay(5.4e-3);
		}
    else  if ( (d3/2) > (5.4e-3 - pwC4 - 5.75e-6) )
		{
		delay(14.0e-3 - pwN - 4.6e-6 - 5.4e-3 + (d3/2));
		xmtroff();
		obsprgoff();
								/*  5.75us   */
		decshaped_pulse("offC4", pwC4, zero, 0.0, 0.0);
		delay(5.4e-3 - pwC4 - 5.75e-6);
		}

    else
		{
		delay(14.0e-3 - pwN - 4.6e-6 - 5.4e-3 + (d3/2));
		xmtroff();
		obsprgoff();
		delay(5.4e-3 - (d3/2) - pwC4 - 5.75e-6);
								/*  5.75us   */
		decshaped_pulse("offC4", pwC4, zero, 0.0, 0.0);
		delay(d3/2);
		}

	rlpower(tpwr, TODEV);					/*  2.3us    */
								/*  point j  */
	sim3pulse(pw, 0.0, pwN, zero, zero, t10, 0.0, 0.0);

	dec2phase(zero);
	delay(2.4e-3 - 1.5*pwN);

	sim3pulse(2.0*pw, 0.0, 2.0*pwN, zero, zero, zero, 0.0, 0.0);

	rlpower(dpwr2, DO2DEV);					/*  2.3us    */
	delay(2.4e-3 - 2.3e-6 - pwN);
		

status(C);
	setreceiver(t11);

}		 
