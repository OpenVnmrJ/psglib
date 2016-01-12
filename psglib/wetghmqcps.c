/* ghmqcps - a phase-sensitive PFG HMQC
   PS version: uses three gradients - all set separately

	j = 1JXH in Hz (140 typical for 1H-13C)
	pwxpwr = (NOT PWXLVL!!) decoupler pulse power level
	pwx = decoupler pulsed pw90 
        gzlvl1 = gradient amplitude (-32768 to +32768)
        gt1 = gradient time (duration) in seconds (0.001)
        gzlvl2 = gradient amplitude (-32768 to +32768)
        gt2 = gradient time (duration) in seconds (0.001)
        gzlvl3 = gradient amplitude (-32768 to +32768)
        gt3 = gradient time (duration) in seconds (0.001)
        grise = gradient rise and fall time (in seconds; 0.00001)
        gstab = optional delay for stability (in seconds)
        nt - works with nt=1 (nt=2 improves data)
        phase - use phase=1,2 to select N,P-type selection to be sorted later
                (use phase=1 to generate an absolute value dataset)

   gzlvl1, gzlvl2, and gzlvl3 and their times (gt1,gt2,gt3) may eventually be fixed
   in their relationship (i.e.2:2:-1, etc) and driven off of one parameter

  Gradient Recommendations:
        for 13C try:
		gzlvl1=20000
                gt1=0.002
                gzlvl2=20000
                gt2=0.002
                gzlvl3=10050
                gt3=0.002

        for 15N try: 
		gzlvl1=
                gt1=0.00
                gzlvl2=
                gt2=0.00
                gzlvl3=00
                gt3=0.00

  PROCESSING:

        process phase=1 (N-type) data with wft2d(1,0,0,1)
        process phase=2 (P-type) data with wft2d(1,0,0,-1)
                   the ('t2dc') argument to wft2d may be useful

        process phase sensitive data (phase = 1,2) with:
                wft2d(1,0,0,1,0,1,1,0) (wft2dnp)

	PAK 920421 - standardized ghmqc.c
	PAK 920430 - ghmqc.c, supports hmbc
	PAK 921209 - ghmqcps.c (tested for 13C organics)
	PAK 930128 - modified
	PAK 930326 - user needs to fiddle with D2_FUDGEFACTOR to empirically make lp1=0
	PAK 940308 - played with d2 timing
	PAK 950522 - added BJ-style phase-cycling,placed grads inside,d2 fiddled
		(tested OK in 3mm; nt=1,2,8,and nt=4 HMBC)
		(if gradients are ++-; used wft2d(10100-101)
	PAK 950523 - made wet version
		(tested OK on LC: nt=8 & 32 (HMBC nt=32 & 8) wft2d(10100-101))
	P.A.Keifer 950920 - updated wet
        P.A.Keifer 960116 - added tpwrf control to wet4
*/


#include <standard.h>
#define D2_FUDGEFACTOR 1.0e-06      /* to empirically make lp1=0 */
static double d2_init = 0.0;

static int ph1[8] = {0,0,0,0,0,0,0,0};
static int ph2[8] = {0,0,0,0,0,0,0,0};
static int ph3[8] = {0,0,2,2,0,0,2,2};
static int ph4[8] = {0,0,0,0,0,0,0,0};
static int ph5[8] = {0,0,0,0,0,0,0,0};
static int ph6[8] = {0,0,0,0,0,0,0,0};
static int ph7[8] = {0,0,0,0,0,0,0,0};
static int ph8[8] = {0,0,0,0,0,0,0,0};
static int ph9[8] = {0,2,0,2,0,2,0,2};
static int ph10[8] = {0,2,2,0,0,2,2,0};
 

pulsesequence()
{
  double j,pwxpwr,gzlvl1,gt1,gzlvl2,gt2,gzlvl3,gt3,grise,gstab,phase,taumb;
  int  icosel,t1_counter;
  char mbond[MAXSTR];

  sw1 = getval("sw1");
  j = getval("j");
  pwxpwr = getval("pwxpwr");
  pwx = getval("pwx");
  gzlvl1 = getval("gzlvl1");
  gt1 = getval("gt1");
  gzlvl2 = getval("gzlvl2");
  gt2 = getval("gt2");
  gzlvl3 = getval("gzlvl3");
  gt3 = getval("gt3");
  grise = getval("grise");
  gstab = getval("gstab");
  phase = getval("phase");
  getstr("mbond", mbond);
  if (mbond[0] == 'y')
      taumb = getval("taumb");

  tau=1/(2.0*j);

  if(tau < (gt3+grise)) 
  {
    text_error("tau must be greater than gt3+grise\n");
    abort(1);
  }

  settable(t1,4,ph1);
  settable(t2,8,ph2);
  settable(t3,4,ph3);
  settable(t4,4,ph4);
  settable(t5,4,ph5);
  settable(t6,4,ph6);
  settable(t7,4,ph7);
  settable(t8,4,ph8);
  settable(t9,4,ph9);
  settable(t10,4,ph10);

/* Gradient incrementation for hypercomplex data */

   if (phase == 2)     /* Hypercomplex in t1 */
   {
      icosel = -1; /* change sign of gradient */
   }
   else icosel = 1;

/* calculate modification to phases based on current t1 values
   to achieve States-TPPI acquisition */
 
   if(ix == 1)
      d2_init = d2;
 
      t1_counter = (int) ( (d2-d2_init)*sw1 + 0.5);
       
      if(t1_counter %2)  
      {
        tsadd(t9,2,4);
        tsadd(t10,2,4);
      }
 

  status(A);
     decpower(pwxpwr);
     delay(d1);
   if (getflag("wet")) wet4(zero,one);
     decpower(pwxpwr);
     rcvroff();

  status(B);
     rgpulse(pw,t1,rof1,rof2);
     delay(tau);
      if (mbond[0] == 'y')                      /* one-bond J(CH)-filter */
      {
         decrgpulse(pwx, t2, rof1, 0.0);
         delay(taumb - tau - rof1 - pwx);
      }

   decrgpulse(pwx,t3,rof1,rof2);
     decphase(t5);
     delay(gt1 + (grise*2.0) + (2.0*GRADIENT_DELAY));
     decrgpulse(pwx*2.0,t3,rof1,rof2);
     txphase(t6);
     rgradient('z',gzlvl1);
     delay(gt1+grise);
     rgradient('z',0.0);
     delay(grise);
   if (d2 > 0);   
     delay(d2/2 + (2.0*pwx/3.14159) - pw + D2_FUDGEFACTOR);

 rgpulse(pw*2.0,t6,rof1,rof2);

     decphase(t9);
     rgradient('z',gzlvl2);
     delay(gt2+grise);
     rgradient('z',0.0);
     delay(grise);
   if (d2 > 0);   
     delay(d2/2 + (2.0*pwx/3.14159) - pw + D2_FUDGEFACTOR);
     decrgpulse(pwx*2.0,t9,rof1,rof2);
     decphase(t9);
     delay(gt2 + (grise*2.0) + (2.0*GRADIENT_DELAY));
     decrgpulse(pwx,t9,rof1,rof2);

     rgradient('z',(double)icosel*gzlvl3);
     delay(gt3+grise);
     rgradient('z',0.0);
     decpower(dpwr);
     setreceiver(t10); 
     rcvron();
     delay(tau-(gt3+grise));
/*     delay(gstab); (Needs to be zero to phase f2 automatically) */
 
  status(C);
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

