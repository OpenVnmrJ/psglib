 /*   dipshft - cross polarization dipshift

 sequence uses no refocussing on the decoupler
 dipolar evolution is via mrev-8

  evan williams, palo alto, october 1989.

Parameter definitions
	dm = 'nny'
	dmm = 'c'
	level1 = decoupler power for cross-polarization (in dhp units)
        level1f = decoupler power for cross-polarization
	level2= decoupler power for dipolar decoupling (in dhp units)
        level2f = decoupler power for dipolar decoupling
	p2 = contact time
	srate= spinning speed
	xpol = 'y' cp expt, = 'n' simple s2pul
	pp = 1H 90 at dipolr power
	pw = 1H 90 at crossp power

 */
#include "standard.h"

mrev8sw(tau,dpuls) double tau,dpuls;
/* semi windowless mrev-8,
     burum et al, jmr, 44, p173, (1981). */

{   


decphase(three);
decon();
delay(dpuls);
decphase(zero);
delay(dpuls);
decoff();
delay(tau);
decphase(zero);
decon();
delay(dpuls);
decphase(one);
delay(dpuls);
decoff();
delay(tau);
decphase(three);
decon();
delay(dpuls);
decphase(two);
delay(dpuls);
decoff();
delay(tau);
decphase(two);
decon();
delay(dpuls);
decphase(one);
delay(dpuls);
decoff();
/* delay(3.0*tau + 8.0*dpuls); */
} 






pulsesequence()
{
  double pp,tau,level1,level1f,level2,level2f,p2,srate,dw,dec;
  char xpol[MAXSTR];
   
    level1=getval("level1");
    level1f=getval("level1f");
    level2=getval("level2");
    level2f=getval("level2f");
    p2=getval("p2");
    srate=getval("srate");
    getstr("xpol",xpol);
    pp=getval("pp");
    tau=getval("tau");
    initval(level1f,v8);
    initval(level2f,v10);
    if (xpol[0] == 'n') 

    {   /* gated decoupling */
    
    /* equilibrium period */
        status(A);
        hsdelay(d1);
    
    /* tau delay */
        status(B);
        pulse(p1,zero);
        hsdelay(d2);
    
    /* observe period */
        status(C);
        power(v9,DODEV);
        pwrf(v10,DODEV);
        obspulse();
    
    } /* end gated decoupling section */
    else /* cross polarization section */
    {   
    
    /* test for any bad parameters */
    
        if (srate < 100.0 ) 
        {   
        fprintf(stdout,"srate must be specified\n");
            abort();
        } 
    
        dw = 8.0*pp + 4.0*tau;
        if (ix == 1) fprintf(stdout,"sw = %6.1f \n",1.0/dw);
        dec = 2.0*pw+p2+np/2.0*dw;
        if (tau < 0.0 ) 
        {   
        fprintf(stdout,"sw too large for mrev!  decrease pp.\n");
            abort();
        } 
    /* check duty cycle for decoupler */
    
        if (dec/(dec+d1) > 0.2 ) 
        {   
        fprintf(stdout,"decoupler dutycycle too high!\n");
            abort();
        } 
    /* check that dm[0]<>y */
    
        if (dm[0] == 'y') 
        {   
        fprintf(stdout,"do not use dm=y..\n");
            abort();
        } 

    /* set up the power levels */

	initval(level1,v7);
	initval(level2,v9);
	initval( (double) (ix-1),v5);

    /* do the pulse sequence */
    
        status(A);
        power(v7,DODEV);	/* set crospolarization power on dec */
        delay(d1+1.0e-6);       /* make sure there is some delay */
	rcvroff();
	delay(rof1);
        decpulse(pw,one);       /* 90(h) */
        decphase(zero);
        status(C);       /* spinlock */
    
    /* contact */
    
        rgpulse(p2,zero,0.0,0.0);
    /* observe section */
    
        if (level2 != level1) power(v9,DODEV);	/* decouple power  level */

        if (ix >= 2)	/*mrev after 1st increment*/
        {
         status(A); 
         if (ix > 2) starthardloop(v5);	/* loop if more than 1 time thru */
         delay(tau);
         mrev8sw(tau,pp);
         if (ix > 2) endhardloop();
        }
        status(C);
        delay(rof2);       /* receiver deadtime */
        rcvron();
    } 
} 
