/* VARIAN VXR 400S COMPOSITE 180 + OBSERVE PULSE SEQUENCE: SANC1 
   Offset alternating Sequence as per Lee et. al. "Multipole Theory of
   Composite Pulses II:Time Reversal and Offset Alternating Sequences
   With Experimental Confirmation", Journal of Magnetic Resonance, 92,455, 1991 */

/*June 26, 1989, Fernando Commodari*/
/*modified September 20,1989;Revised September 20, 1991, F.C.
        p1    p1    p1                  pw
       ____  ____  ____                ____
  d1   !  !  !  !  !  !            d2  !  !  ACQ
       !  !  !  !  !  !                !  !
------[!f1!--!f2!--!f1!--......]n------!  !-------- 

n: the number of composite pulses to be executed,
   if n=0,a simple observe pulse(pw=p1/2) is executed
   if n=1,a single comp pulse with offset f1 is executed followed by an obs pw
   if n=even, n pulses are executed with f1&f2 being looped (n/2) times
      followed by the pi/2 observe pulse
   if n=odd not 1,an extra comp pulse with offset f1 is added to  the 
      (n-1) even case followed by the observe pulse  
p1:the length of the individual composite pulses;***set p1=2*pw/n***
pw:the observe pulse
d1:relaxation delay
d2:optional evolution delay
to:desired offset from resonance frequency,tof
f1:odd pulse offset (pulse 1,3,5...) calculated:f1=(to+tof)>0,if tof=0
f2:even pulse offset(pulse 2,4,6,8....)calculated:f2=(tof-to)<0,if tof=0*/


#include <standard.h>

pulsesequence()
{
double f1,f2,n,to;
n=getval("n");
to=getval("to");
initval(n,v10);
initval(2.0,v5);
f1=(to+tof);
f2=(tof-to);
mod2(v10,v2);
decr(v10);
assign(v10,v3);
divn(v10,v5,v4);

/*Preparation Period*/

  status(A);
  delay(d1);

/*Evolution Period*/

  ifzero(v10);/* if n=0,execute simple observe (pi/2) pulse */
   status(B);
   offset(tof,TODEV);
   delay(d2);
   status(C);
   rgpulse(pw,oph,rof1,rof2);
  elsenz(v10);
  ifzero(v3);/* if n=1, execute one comp pulse with offset f1 */
   status(B);
   offset(f1,TODEV);rgpulse(p1,oph,0.0,0.0);
  elsenz(v3);/* for even # of comp pulses,n, loop f1&f2 (n/2) times */
   status(B);
   starthardloop(v4);
   offset(f1,TODEV);rgpulse(p1,oph,0.0,0.0); 
   offset(f2,TODEV);rgpulse(p1,oph,0.0,0.0);
   endhardloop();
  ifzero(v2);
  elsenz(v2);/* even/odd flag,if v2 not 0,n=odd not 1,add extra f1 to (n-1) even  case*/ 
   offset(f1,TODEV);rgpulse(p1,oph,0.0,0.0);
  endif(v2);
  endif(v3);
  endif(v10);
offset(tof,TODEV);
delay(d2);
status(C);
pulse(pw,oph);
}



 
