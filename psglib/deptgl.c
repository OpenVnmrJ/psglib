#ifndef LINT
static char SCCSid[] = "@(#)deptgl.c 1.8 1/7/89    Copr 1988 Varian Assoc.";  
#endif

/* deptgl - for spectral editing and polarization transfer experiments */ 

/* parameters:
      pw	observe 90
      pp	decoupler 90
      pplvl     power level for decoupler proton pulses (for 500);
      d2,d3,d4	j evolution delays, operator enterable if jmax and jmin
 		are 0, or calculated if jmin<>0 and jmax<>0
      jmax,jmin	max and min j values, used to calculate d2,d3,d4 to
		minimize J cross talk
      degree	editing flip angle
		for editing set to 38,90,90,142
      satdly    optional saturation delay (ca. 0.05 sec)
      nt	set to 64 transients minimum      */

/* reference: 
      O.W.Sorensen et al., J.Magn.Reson. 55, 347 (1983). */


#include <standard.h>


pulsesequence()
{
   double ppp,d3,d4,d5,pp,ps,jmin,jmax,degree,satdly,pplvl;
   
/* get new variables */
   pp = getval("pp");
   d3 = getval("d3");
   d4 = getval("d4");
   jmin = getval("jmin");
   jmax = getval("jmax");
   degree = getval("degree");
   satdly = getval("satdly");
   dfrq = getval("dfrq");
   if (newdecamp)
      {
       pplvl = getval("pplvl");
       initval(pplvl,v14);
       initval(dhp,v13);
       if (rof2 == 0.0)
          {
           rof2 = 1.0e-5;
          }
       ps = 0.4e-6;
      }
   else
      {
       ps = 5.0e-6;
      }
   if (rof1 < 1.0e-5)
      {
       rof1 = 1.0e-5;
      }

/*setup calculations*/
   if ((((jmax == 0.0) || (jmin == 0.0)) &&
          ((d2 == 0.0) || (d3 == 0.0) || (d4 == 0.0))) ||
         (jmin > jmax))
      {
         fprintf(stdout,
            "wrong parameters: jmin and jmax reset to 125 and 200 Hz.\n");
 	jmin = 125.0; jmax = 200.0;
      }
   if ((jmin != 0.0) && (jmax != 0.0))
      {
         ppp = (jmax-jmin)*0.07;
	 d2 = 1.0/(2.0*(jmin+ppp));
	 d3 = 1.0/(2.0*(jmax-ppp));
	 d4 = 1.0/(jmax+jmin);
         fprintf(stdout,
            "degree = %5.1f:  d2 = %f, d3 = %f, d4 = %f\n",degree,d2,d3,d4);
      }
   if (d4 > d2)
      {
 	d5 = d4; d4 = d2; d2 = d5;
      }
   if (degree > 90.5)
      {
         degree = 180.0-degree; assign(two,v1);
      }
   else
      {
       assign(zero,v1);
      }
   ppp = (pp*degree)/90.0;


/* pulse sequence code */

      if ((dm[0] == 'y') || (dm[1] == 'y')) 
         {   
            fprintf(stdout,"decoupler must be set as dm='nny'\n");
            abort(1);
         } 
    if ((dmm[0] != 'c') || (dmm[1] != 'c'))
      {
       fprintf(stdout,"Decoupler must be set as dmm='cce' or dmm='ccs'.\n");
       abort(1);
      }

/* ACTUAL PULSESEQUENCE BEGINS */
   status(A);
       if (dfrq > 495)
         {
          delay(0.2e-6); power(v14,DODEV);
         }
       else
         {
          declvlon();   /* sets dhp=255 level for non-VXR500 RF architecture */ 
         }
       hsdelay(d1);
       if (satdly > 0.0 )
         {
          status(B);
          pulse(pw,zero);
          hsdelay(satdly);       /* destroy xy magnetization */
         }
     
/* setup phase cycle */
   sub(three,oph,v7);	/* 3210 */
   incr(v7);            /* 0321 */
   mod2(ct,v6);         /* 0101 */
   dbl(v6,v6);		/* 0202 */
   hlv(ct,v9);		/* 0011 2233 */
   hlv(v9,v9);		/* 0000 1111 2222 3333 */
   mod2(v9,v4);         /* 0000 1111 0000 1111 */
   dbl(v4,v4);		/* 0000 2222 0000 2222 */
   hlv(v9,v9);          /* ct/8 */
   mod2(v9,v5);         /* 0000 0000 1111 1111 */
   dbl(v5,v5);		/* 0000 0000 2222 2222 */
   add(v4,v7,v7);	/* 0321 2103 0321 2103 */
   add(v5,v7,v7);	/* 0321 2103 2103 0321 */
   add(v1,v5,v3);	/* 0000 0000 2222 2222 or */
   			/* 2222 2222 0000 0000 */
   add(v6,v5,v5);	/* 0202 0202 2020 2020 */
   mod4(v9,v8);		/* ct/8 */
   hlv(v9,v2);		/* ct/16 */
   mod2(v2,v9);         /* 0000 0000 0000 0000 1111 1111 1111 1111 */
   dbl(v9,v9);		/* 0000 0000 0000 0000 2222 2222 2222 2222 */
   add(v4,v9,v9);	/* 0000 2222 0000 2222 2222 0000 2222 0000 */
   incr(v9);		/* 1111 3333 1111 3333 3333 1111 3333 1111 */
   add(v6,v9,v9);	/* 1313 3131 1313 3131 3131 1313 3131 1313 */

   status(B);
      decrgpulse(pp,zero,rof1,rof2);
      if ((0.5*(d2-d4))<(pp+2e-7))
         {
          delay(d4-pp-rof1-rof2);
          simpulse(pw,2.0*pp,v7,v2,rof1,rof2);
         }
      else if ((0.5*(d2-d4)-pp-rof1-rof2)<2e-7)
         {
          delay(d4-pw-rof1-rof2);
          rgpulse(pw,v7,rof1,0.0);
          decrgpulse(2.0*pp,v2,0.5*(d2-d4)-pp,rof2);
         }
      else
         {
          delay(d4-pw-rof1-rof2);
          rgpulse(pw,v7,rof1,rof2);
          delay(0.5*(d2-d4)-pp-rof1-rof2);
          decrgpulse(2.0*pp,v2,rof1,rof2);
         }
      delay(0.5*(d2+d4)-pp-rof1-rof2);
      decrgpulse(pp,v3,rof1,rof2);
      delay(d3/2.0-pp-rof1-rof2);
      simpulse(2.0*pw,2.0*pp,v8,oph,rof1,rof2);
      delay(d3/2.0-pp-rof1-rof2);
      decrgpulse(pp,v5,rof1,0.0);
      decrgpulse(ppp,v9,ps,rof2);
    
      if (dm[2] == 'n')
         {
            delay(d2-ppp-2.0*pp-rof1-ps-rof2);
	    ifzero(v6);
               decrgpulse(pp,zero,rof1,0.0);
               decrgpulse(2.0*pp,one,ps,0.0);
               decrgpulse(pp,zero,ps,rof2);
	    elsenz(v6);
	       delay(4.0*pp+2.0*ps+rof1+rof2);
	    endif(v6);
	 }
      else
	 {
            delay(d2-rof2);
  	 }
      if (dfrq > 495)
         {
          delay(0.2e-6); power(v13,DODEV);
         }
      else
         {
          declvloff();
         }
      decphase(zero);
   status(C);
} 
