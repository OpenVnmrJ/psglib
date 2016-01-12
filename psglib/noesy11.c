/* noesy11 -  2D cross relaxation experiment using 11echo "read" pulse
 
       mix = mixing time.
     phase = 1,2: gives HYPERCOMPLEX phase-sensitive experiment;
        p1 = correction of second pulse width for radiation damping
        p2 = correction of second pulse width for radiation damping
       tau = 1-1 pulse separation(typ 50-100us: try 80us)
      tau2 = times flanking 1-1 echo part  (<<1ms) (try 100us)
       hst = 40-60ms or half the mix time.Necessary to eliminate x-y H2O
             component prior to 11echo read pulse. 
        hs = nnyn or ynyn
    
     use "calfa" to set alfa for lp=0 to obtain flat baseline

                                Sklenar and Bax, JMR, 74,469(1987)  */
#include <standard.h>
pulsesequence()
{
   double mix,tau,p2,tau2;
   p2=getval("p2"); tau=getval("tau"); mix=getval("mix");
   tau2 = getval("tau2"); loadtable("noesy11");
   ttadd(t1,t8,4); ttadd(t2,t8,4); ttadd(t3,t8,4);    /* add cyclops */
   ttadd(t4,t8,4); ttadd(t5,t8,4); ttadd(t6,t8,4); ttadd(t7,t8,4);

   getelem(t1,ct,v1);
   if (getval("phase") > 1.5) incr(v1);      /* hypercomplex phase shift */

   status(A);
     hsdelay(d1);
   status(B);
      rgpulse(pw, v1, rof1, 1.0e-6);
      delay(d2);
      rgpulse(pw, t2, rof1, 1.0e-6);
   status(C);
      rcvroff();
      hsdelay(mix);
   status(D);
      rgpulse(pw,t3,rof1,rof1);
      delay(tau-pw-rof1);
      rgpulse(pw-p1,t4,rof1,rof1);
      hsdelay(tau2-rof1-pw);
      rgpulse(pw,t5,rof1,rof1);
      delay(2.0*tau-pw-rof1);
      rgpulse(pw-p2,t6,rof1,rof1);
      hsdelay(tau2-rof1-pw/2.0);
      setreceiver(t7);
      rcvron();
                         /* noesy11 tablib entry
                          t1 = (0)32
                          t2 = (0 0 0 0 2 2 2 2)4 
                          t3 = (0)32 
                          t4 = (2)32
                          t5 = (0 1 2 3)8
                          t6 = (2 3 0 1)8
                          t7 = (0 2 0 2 2 0 2 0)4
                          t8 = [0 1 2 3]8    */
}
