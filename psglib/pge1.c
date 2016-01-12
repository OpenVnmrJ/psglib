/*pgeramp1.c - diffusion experiment with trapezoidal ramped gradients
               from DIFFUSION.tar.Z 
               D. Rice - Modified 11/27/06 to add psg_abort */

#include <standard.h>
#include <math.h>

/*function to create a trapezoidal gradient pulse*/ 

rampgrad(on_val,off_val,length,ramp)
double on_val,off_val,length,ramp;
{
   int       steps,i;

   double    initval,iramp,incr_val;

   if (ramp > length) {
       fprintf(stdout,"ramp parameter larger than gradient delay\n");
       psg_abort(1);
   }

   steps = (int) (ramp /(1.0e-6 + GRADIENT_DELAY));
   if (steps > 100)  steps = 100;  

   incr_val=0.0;
   initval=0.0;
   iramp=0.0;
   if (steps > 0) {
      incr_val = (on_val - off_val) / steps;
      initval = off_val + incr_val;
      iramp = (ramp / (double) steps) - GRADIENT_DELAY;
   }
   
   if (length >= 2e-7) {     
      if (steps > 0) {
         for (i = 0; i < steps; i++) {
            rgradient('Z', initval);
            delay(iramp);
            initval = initval + incr_val;
         }
      }
      rgradient('Z',on_val);
      delay(length - ramp);
      initval = on_val - incr_val; 
      if (steps > 0) {   
         for (i = 0; i < steps; i++) {
            rgradient('Z', initval);
            delay(iramp);
            initval = initval - incr_val;
         }
      }
      rgradient('Z',off_val);
   }
} 

static int table1[4] = {0, 1, 2, 3};
static int table2[4] = {0, 1, 0, 1};

pulsesequence()
{

/* declare new variables*/

     double   dutycycle,p2,g1,g2,d0,d5,dac_cw,at,grad_cw,grad_p1,dac_p1,
              time,tramp;

     char     stim[MAXSTR];    

/*get new variables*/

     p2=getval("p2");
     g1=getval("g1");
     g2=getval("g2");
     d0=getval("d0");
     d5=getval("d5");
     dac_cw=getval("dac_cw");
     grad_cw=getval("grad_cw");
     grad_p1=getval("grad_p1");
     dac_p1=getval("dac_p1");
     at=getval("at");
     tramp=getval("tramp");
     getstr("stim",stim);

/*check the gradient duty cycle*/

     time=g1+g2+d0+d1+d2+d3+d4+d5+p1+p2+pw+at;
     dutycycle=sqrt(((grad_cw*grad_cw*time)+(grad_p1*grad_p1*(g1+g2)))/time);
     
     if (dutycycle>100) {
        fprintf(stdout,"ABORT: Average Gauss/cm = %g. Gauss/cm must be < 100\n",
                  dutycycle);
        psg_abort(1);
     }

     if ((dac_cw > 0) && (ix == 1)) {
        if (dac_cw > 100) {
           fprintf(stdout,"ABORT: CW Gauss/cm is too large\n");
           psg_abort(1);
        }     
        fprintf(stdout,"WARNING: CW Gradient is ON at %g guass/cm\n" ,grad_cw);
     }

/*Set Phase Tables*/

     settable(t1,4,table1);
     settable(t2,4,table2);
     setreceiver(t1);

/*begin pulse sequence*/

     status(A);
     txphase(t1); 
     rgradient('Z', dac_cw);
     delay(d1);
     sp1on(); delay(2.0e-6); sp1off(); delay(2.0e-6);
     rcvroff();
     delay(rof1);    
     rgpulse(p1,t1,0.0,0.0);
     status(C);
     delay(d0);
     dps_off();
     rampgrad(dac_p1,dac_cw,g1,tramp);
     dps_on();
     dps_show("zgradpulse",dac_p1,g1-tramp); 

     if(stim[0] == 'y') {
        txphase(t1);
        delay(d2);
        rgpulse(p2,t1,0.0,0.0);
        delay(d5);
        rgpulse(pw,t1,0.0,0.0);
     }
     else {
        txphase(t2); 
        delay(d2);                   
        rgpulse(pw,t2,0.0,0.0);
     }
     delay(d3);
     dps_off();
     rampgrad(dac_p1,dac_cw,g2,tramp);
     dps_on();
     dps_show("zgradpulse",dac_p1,g2-tramp); 
     delay(d4);

/*begin acquisition*/

     obsblank();
     delay(rof2);
     rcvron();
}
