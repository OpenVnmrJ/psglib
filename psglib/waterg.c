/* 1D watergate version 2  RC 31 Oct 94 

	psel	1H softpulse	(5350)
	psellvl	pwr for psel	(-1)
	phaseinc  0.5 degree phaseshift for psel = (1)

 Calibrate a 5 to 6 ms 1H softpulse. After you have decided upon
 tpwr and pw(90) refine with sequence.  Set psel to ~ 10 microsec
 accuracy.
 MILD DPWR and DPWR2 PARAMETER CHECKING IS DONE THOUGH IT IS
 YOUR RESPONSIBILITY TO NOT MELT YOUR PROBE. 
 Set dm and/or dm2 = 'nny' to decouple 13C and 15N with labeled
 samples.

R. Crouch Burroughs Wellcome Co.
crouch@bwco.com				*/


# include <standard.h>

static int	phs1[4] = {0,2,1,3},
		phs2[4] = {2,0,3,1};

pulsesequence()
{
 double psel,psellvl,gt1,glvl1,grise,phaseinc;

 psel = getval("psel");
 psellvl = getval("psellvl");
 gt1 = getval("gt1");
 glvl1 = getval("glvl1");
 grise = getval("grise");
 phaseinc = getval("phaseinc");

/* Check for correct DM settings */
   if ((dm[A] == 'y') || (dm2[A] == 'y'))
   {
      printf("DM/DM2 must be set to either 'nny' or 'nnn'.\n");
      abort(1);
   }

/* Check dpwr & dpwr2 */
   if ((dpwr > 52) || (dpwr2 > 50))
   {
	printf("DPWR and/or DPWR2 out of range. Don't melt the probe.\n");
	abort(1);
   }

 
 if (phaseinc < 0.0) 
   phaseinc = 720+phaseinc;
 stepsize(0.5, TODEV);
 initval(phaseinc, v1);

 settable(t1,4,phs1);
 settable(t2,4,phs2);

 status(A);
 delay(d1);
 rcvroff();
 status(B);
 pulse(pw,t1); delay(rof1);
 delay(grise);
 zgradpulse(glvl1, gt1);
 rlpower(psellvl, TODEV);
 xmtrphase(v1);
 delay(grise);
 txphase(t2);
 pulse(psel,t2);
 xmtrphase(zero);
 txphase(t1); rlpower(tpwr, TODEV);
 delay(2e-6);
 pulse(2*pw,t1);
 txphase(t2);
 rlpower(psellvl, TODEV); 
 xmtrphase(v1);
 delay(2e-6);
 pulse(psel,t2);
 delay(grise);
 zgradpulse(glvl1, gt1);
 rlpower(tpwr, TODEV);
 xmtrphase(zero);
 delay(grise);
 setreceiver(t1);
 status(C);
 rcvron(); delay(rof1);
}
