/* xdec2pul.c  Two-pulse sequence in "reverse" configuration 
               Uses H1 decoupler as observse and BB transmitter as decoupler
 

	Parameters:
   	dn,dof 	used to set proton (observe) freq
	pw	90 degree pulse for protons
	dpwr	pulse power for H1 pulse
	tn,tof	used to set X nucleus freq
	tpwr	power level for X nucleus pulse
	xdec    two states, yn gives  X nucleus decoupling 
		during D1 for NOE studies
                ny gives decoupling during acquisition
                n gives coupled spectrum
    modified 6/12/92 Steve Cheatham - removed status & fixed
    dm=ynn problem 	*/

#include <standard.h>
#include <revmode.c>

pulsesequence()
{   
   double 	dx1;
   char		xdec[MAXSTR];
   
   getstr("xdec",xdec);   
   

   dx1=0.88/sw; 	/* set 12% decoupling duty cycle */
   if (((1.0/sw)-dx1) < 4.0e-7)
    {
     printf("xdec2pul: SW too large for decoupling during acquisition.\n");
     abort(1);
    }
  if ((dm[A]=='y') || (dm[B]=='y') || (dm[C]=='y')) /* force dm to nnn */
    {
     printf("Cannot use dm in this pulse sequence, dm must = nnn.\n");
     abort(1);
    }

   if (xdec[A]=='y')  { rcvroff(); xmtron(); }    
   hsdelay(d1);
  xmtroff(); 
  decrgpulse(p1,oph,rof1,rof2);
  hsdelay(d2);
  decrgpulse(pw,oph,rof1,rof2);

  if (xdec[B]=='y')	/* check for X decoupling flag */
  {			/* xdec='y' decoupler on */
   delay(alfa+1.0/(beta*fb));	/* receiver gate delay */
   acquire(2.0,2.0e-7);		/* acquire 1st set of points */
   initval((np/2)-1.0,v9);
   starthardloop(v9);		/* acquire the rest of np points */
     rcvroff();			/* gate receiver off */
     xmtron();			/* gate decoupler on */
     delay((1.0/sw)-dx1 - 2.0e-7); /* time shared decoupling */
     xmtroff();			/* gate decoupler off */
     rcvron();			/* gate receiver on */
     delay(dx1);		
     acquire(2.0,2.0e-7);	/* acquire 1 set of points */
   endhardloop();
  }				/* xdec='n', implicit acquisition on */
} 
