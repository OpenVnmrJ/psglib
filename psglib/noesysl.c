
/* noesysl -  2D cross relaxation experiment using 
              trim pulses for water suppression and phase tables 
               G.Otting et.al., J.Biomolecular NMR,1,209(1991)

               typical values d3=60us  slock1=500us slock2=2000us
               slpower=full to 6-12 db lower(test water     
               suppression at different slpower and use minimum
               value to minimize sample heating.

               (use of full power for trim pulses will usually
                increase sample temperature. allow enough dummy
                pulses to achieve equilibrium. this could be
                several minutes.)

*/

#include <standard.h>

pulsesequence()
{
   double mix,slpower,slock1,slock2; 
   int iphase;
   loadtable("noesysl");

/* LOAD VARIABLES */
   mix = getval("mix"); slpower = getval("slpower");
   slock1 = getval("slock1"); slock2 = getval("slock2");
   iphase = (int) (getval("phase") + 0.5);

/* STATES-HABERKORN */
   getelem(t1,ct,v1);  
   if (iphase == 2)
      incr(v1);
   getelem(t6,ct,oph);

/*HYPERCOMPLEX MODE USES REDFIELD TRICK TO MOVE AXIAL PEAKS TO EDGE */
    initval(2.0*(double)(((int)(d2*getval("sw1")+0.5)%2)),v6);
  add(v1,v6,v1); add(oph,v6,oph);  

/* BEGIN THE ACTUAL PULSE SEQUENCE */
   status(A);
      rlpower(tpwr,TODEV);
      hsdelay(d1);
      rgpulse(pw, v1, rof1, 1.0e-6);
      if (d2>0.0) delay(d2 - rof1 - 1.0e-6 - (4.0*pw/3.14159));
      rgpulse(pw, t2, rof1, 1.0e-6);
      delay(mix);
      rgpulse(pw, t3, rof1, rof1);
      rlpower(slpower,TODEV);
      rgpulse(slock1,t4,rof1,rof1);
      delay(d3);
      rgpulse(slock2,t5,rof1,rof2);
}
/* noesysl table 

t1 = [0 2]8
t2 = (0 2)16
t3 = (0 0 2 2)8 
t4 = [1 3]4
t5 = [0 2]16
t6 = (0 2 2 0)2 (2 0 0 2)2
*/
