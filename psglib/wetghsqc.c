/* ghsqc_da.c 
   heteronuclear single quantum correlation with gradient selection and editing.

    THIS PULSE SEQUENCE DESIGNED TO USE CHANNELS 1 and 2.

    gt1     - first gradient duration (0.002 for 13C)
    gzlvl1  - first gradient power level (20000 for 13C)
    gt2     - second gradient duration (0.0005 for 13C)
    gzlvl2  - second gradient duration (-19800 fro 13C - optimize!)
	      The gradients need to maintain a 4:1 (13C) or 10:1 (15N) ratio:
	      the ratio can be controlled with either amplitude or time.
    pwxpwr  - (NOT pwxlvl!) decoupler power level for hard decoupler pulses
    pwx     - pulse width for hard decoupler pulses
    dpwr    - power level for decoupling
    dmf     - controlled by dpwr
    dmm     - 'ccp'
    dseq    - mpf7 (covers ca. 8 x YB1)
    j       - one-bond heteronuclear coupling constant (140 for 13C; 90 for 15N)
    satmode - (n) a y/n flag for transmitter presaturation
    satdly  - the presaturation delay used if satmode = y
    satpwr  - the presaturation power
    satfrq  - the frequency desired for presaturation
    edit    - (n) a y/n flag for X editing
    fad     - (y) a y/n flag for States-TPPI
    f1180   - (n) a y/n flag for 1/2 dwell starting t1 evolution delay
    nt      - works with nt=1 (nt=2 improves data)
    phase   - use phase=1,2 to select N,P-type selection to be sorted later
              (transform with wft1dnp/wft2dnp [np == (1,0,1,0,0,1,0,-1)]) 

Lit: W. Wilker et al Magn. Reson. Chem. 31,287-292(1993)
	ech jan 95 dast
	P.A.Keifer 950515 - made wetghsqc4.c
		works OK for phase=1,2 (nt=4 OK; nt=1 has FADed junk)
		but lp1 nonzero (d2 goofed up)
	P.A.Keifer 950523 - pwxpwr,d2(removed rof1s),
			presetphases,removed TODEVs and DODEVs
		(tested OK; with lp1 almost zero)
		(process with wft2d(1010010-1); OK with nt=(1),2,4,and 8)
	P.A.Keifer 950920 - updated wet
        P.A.Keifer 960116 - added tpwrf control to wet4
*/

#include <standard.h>
static double d2_init = 0.0;
static int phi1[1] = {0},
           phi2[1] = {0},
           phi3[1] = {0}, 
           phi4[4] = {1,1,3,3},
           phi5[2] = {0,2},
           phi6[1] = {0},
           phi7[2] = {0,2},
           phi8[2] = {0,2},
           phi9[8] = {0,0,0,0,2,2,2,2},
           phi10[16] = {0,0,0,0,0,0,0,0,2,2,2,2,2,2,2,2},
           phi11[1] = {0},
           phi12[1] = {0},
           rec_1[16] = {0,2,2,0,2,0,0,2,2,0,0,2,0,2,2,0};

pulsesequence()
{
/* DECLARE VARIABLES */
 char     f1180[MAXSTR],edit[MAXSTR],fad[MAXSTR];
 int	  phase,satmove,icosel;
 int      t1_counter;  
 double   gt1 = getval("gt1"),
          gt2 = getval("gt2"),
          gzlvl1 = getval("gzlvl1"),
          gzlvl2 = getval("gzlvl2"),
          pwxpwr = getval("pwxpwr"),
          j = getval("j"),
          tau1,tauxh;

  getstr("f1180",f1180);
  getstr("satmode",satmode);
  getstr("edit",edit);
  getstr("fad",fad);

  phase = (int) (getval("phase") + 0.5);
  satmove = ( fabs(tof - satfrq) >= 0.1 );

/* LOAD VARIABLES */
  settable(t1, 1, phi1);
  settable(t2, 1, phi2);
  settable(t3, 1, phi3);
  settable(t4, 4, phi4);
  settable(t5, 2, phi5);
  settable(t6, 1, phi6);
  settable(t7, 2, phi7);
  settable(t8, 2, phi8);
  settable(t9, 8, phi9);
  settable(t10, 16, phi10);
  settable(t11, 1, phi11);
  settable(t12, 1, phi12);
  settable(t13, 16, rec_1);

/* INITIALIZE VARIABLES */
   if (j != 0.0)
      tauxh = 1/(4*(j));
   else
      tauxh = 1.80e-3;
 
/* PHASE INCREMENTATION FOR HYPERCOMPLEX DATA */
   if (phase == 2) 
   {
     icosel = 1;
   }
   else icosel = -1; 

/* CALCULATE MODIFICATION TO PHASES BASED ON CURRENT T1 VALUES
   TO ACHIEVE STATES-TPPI ACQUISITION */
   if(ix == 1)
     d2_init = d2;

   if (fad[A] == 'y')
   {
   t1_counter = (int) ( (d2-d2_init)*sw1 + 0.5);
   
   if (t1_counter %2)
     {
     tsadd(t3,2,4);
     tsadd(t5,2,4);
     tsadd(t13,2,4);
     }
   }

/* SET UP FOR (-90,180) PHASE CORRECTIONS IF f1180='y' */
   tau1 = d2;
   if(f1180[A] == 'y')  tau1 += ( 1.0/(2.0*sw1) );
   tau1 = tau1/2.0;
   
/* BEGIN ACTUAL PULSE SEQUENCE */
 status(A);
   rcvroff();
   decpower(pwxpwr);               /* Set decoupler power to pwxpwr */
   if(satmode[A] == 'y')
   {
   if (satmove)
     obsoffset(satfrq);
     delay(d1-satdly); 
     txphase(zero);
     obspower(satpwr);             /* Set power for presaturation  */
     rgpulse(satdly,zero,rof1,rof2);    /* Presaturation using transmitter */
     obspower(tpwr);               /* Set power for hard pulses  */
     if (satmove)
     obsoffset(tof);
   }
   else
   {
   obspower(tpwr);                 /* Set power for hard pulses  */
   delay(d1);
if (getflag("wet")) wet4(zero,one);
   decpower(pwxpwr);               /* Set decoupler power to pwxpwr */
   }

 status(B);
   rcvroff(); 
   rgpulse(pw,t1,1.0e-6,0.0);    
   txphase(t2); decphase(t3);
   delay(tauxh);                        /* delay=1/4J(XH) */
   simpulse(2*pw,2*pwx,t2,t3,0.0,0.0);
   txphase(t4); decphase(t5);
   delay(tauxh);
   simpulse(pw,pwx,t4,t5,2.0e-6,0.0);	
   if (tau1 > 2.0*pwx/3.1416)
	 delay(tau1 - 2.0*pwx/3.1416 - pw); /* delay=t2/2 */ 
   rgpulse(2*pw,t6,rof1,rof1);
   if (tau1 > 2.0*pwx/3.1416)
	 delay(tau1 - 2.0*pwx/3.1416 - pw);
   if (edit[0] == 'y')
    {
    txphase(t8); decphase(t7);
    delay(2.0*tauxh - 1.0e-4 - gt1 - rof1 - 2.0*GRADIENT_DELAY);
    zgradpulse(icosel*gzlvl1,gt1);
    delay(1.0e-4);
    simpulse(2.0*pw,2.0*pwx,t8,t7,rof1,rof1);
    txphase(t9); decphase(t10);
    delay(2.0*tauxh - rof1);
    }
   else
    {
    decphase(t7);
    zgradpulse(icosel*gzlvl1,gt1);
    delay(1.0e-4);
    decrgpulse(2.0*pwx,t7,rof1,rof1);
    txphase(t9); decphase(t10);
    delay(gt1 + 1.0e-4 + 2.0*GRADIENT_DELAY);
    }
   simpulse(pw,pwx,t9,t10,0.0,0.0);
   txphase(t11); decphase(t12);
   delay(tauxh);                        /* delay=1/4J (XH)  */
   simpulse(2*pw,2*pwx,t11,t12,0.0,0.0);  
   decpower(dpwr);  /* lower decoupler power for decoupling */
   zgradpulse(gzlvl2,gt2);
   delay(tauxh -  gt2 - 2*GRADIENT_DELAY);
   rcvron();
 status(C);
   setreceiver(t13);
}

/* wet4 - Water Elimination */
wet4(phaseA,phaseB)
  codeint phaseA,phaseB;
 
{
  double finepwr,gzlvlw,gtw,gswet,dmfwet,dpwrwet,dofwet,wetpwr,pwwet,dz;
  int c13wet;
  char   wetshape[MAXSTR];
  c13wet=getflag("c13wet");             /* Water suppression flag        */  
  getstr("wetshape",wetshape);    /* Selective pulse shape (base)  */
  wetpwr=getval("wetpwr");        /* User enters power for 90 deg. */
  pwwet=getval("pwwet");        /* User enters power for 90 deg. */
  dmfwet=getval("dmfwet");
  dpwrwet=getval("dpwrwet");
  dofwet=getval("dofwet");
  dz=getval("dz");
  finepwr=wetpwr-(int)wetpwr;     /* Adjust power to 152 deg. pulse*/
  wetpwr=(double)((int)wetpwr);
  if (finepwr==0.0) {wetpwr=wetpwr+5; finepwr=4095.0; }
  else {wetpwr=wetpwr+6; finepwr=4095.0*(1-((1.0-finepwr)*0.12)); }
  rcvroff();
  if (c13wet)
    {
    setstatus(DECch,FALSE,'w',FALSE,dmfwet);
    decoffset(dofwet);
    decpower(dpwrwet);
    }
  obspower(wetpwr);         /* Set to low power level        */
  gzlvlw=getval("gzlvlw");      /* Z-Gradient level              */
  gtw=getval("gtw");            /* Z-Gradient duration           */
  gswet=getval("gswet");        /* Post-gradient stability delay */
  chess(finepwr*0.5059,wetshape,pwwet,phaseA,20.0e-6,rof2,gzlvlw,gtw,gswet,c13wet);
  chess(finepwr*0.6298,wetshape,pwwet,phaseB,20.0e-6,rof2,gzlvlw/2.0,gtw,gswet,c13wet);
  chess(finepwr*0.4304,wetshape,pwwet,phaseB,20.0e-6,rof2,gzlvlw/4.0,gtw,gswet,c13wet);
  chess(finepwr*1.00,wetshape,pwwet,phaseB,20.0e-6,rof2,gzlvlw/8.0,gtw,gswet,c13wet);
  if (c13wet)
    {
    setstatus(DECch,FALSE,'c',FALSE,dmf);
    decoffset(dof);
    decpower(dpwr);
    }
  obspower(tpwr);               /* Reset to normal power level   */
  obspwrf(tpwrf);
  rcvron();
  delay(dz);
}
/* chess - CHEmical Shift Selective Suppression */
chess(pulsepower,pulseshape,duration,phase,rx1,rx2,gzlvlw,gtw,gswet,c13wet)  double pulsepower,duration,rx1,rx2,gzlvlw,gtw,gswet;
  int c13wet;
  codeint phase;
  char* pulseshape;
{
  obspwrf(pulsepower);
  if (c13wet) decon();
  shaped_pulse(pulseshape,duration,phase,rx1,rx2);
  if (c13wet) decoff();
  zgradpulse(gzlvlw,gtw);

  delay(gswet);
}
 
int getflag(str)
char str[MAXSTR];
{
   char strval[MAXSTR];
 
   getstr(str,strval);
   if ((strval[0]=='y') || (strval[0]=='Y')) return(TRUE);
     else                                    return(FALSE);
}



