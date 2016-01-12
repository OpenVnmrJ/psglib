/* noesyjr -  2D cross relaxation experiment using jr "read" pulse
 
         Normal noesy parameters are used
      Only States method used for phase-sensitive in F1 (no tppi) 

      p1 = correction of second pulse width for radiation damping
      jrdelay = interpulse delay in jr read pulse (try 140us) 
     (use "calfa" to set alfa for lp=0 to obtain flat baseline)

          G.Gray   Palo Alto  Sept 1991
*/
#include <standard.h>
pulsesequence()
{
   double jrdelay,mix;
   char sspul[MAXSTR]; getstr("sspul",sspul);
   mix=getval("mix"); jrdelay=getval("jrdelay");
   loadtable("noesyjr");
   getelem(t1,ct,v1);
   if (getval("phase") > 1.5) incr(v1);      /* hypercomplex phase shift */

   status(A);
     rcvroff();
     if (sspul[A]=='y')
      {
       rgpulse(200*pw,zero,rof1,0.0); rgpulse(200*pw,one,0.0,rof1);
      }
     txphase(v1);
     hsdelay(d1);
   status(B);
      rgpulse(pw, v1, rof1, 1.0e-6);
      txphase(t2);
      if (d2>0.0)
      delay(d2 -rof1 -(4.0*pw)/3.1414);
      rgpulse(pw, t2, rof1, 1.0e-6);
   status(C);
      txphase(t3);
      hsdelay(mix);
   status(D);
      rgpulse(pw,t3,rof1,rof1);
      txphase(t4);
      delay(jrdelay-2.0*rof1-pw);
      rgpulse(pw-p1,t4,rof1,rof1);
      setreceiver(t5);
      rcvron();
}

/*
 noesyjr	

t1 = 0 2 0 2 0 2 0 2 0 2 0 2 0 2 0 2 
t2 = 0 0 0 0 0 0 0 0 2 2 2 2 2 2 2 2 
t3 = 0 2 1 3 2 0 3 1 0 2 1 3 2 0 3 1  
t4 = 2 0 3 1 0 2 1 3 2 0 3 1 0 2 1 3
t5 = 0 0 1 1 2 2 3 3 2 2 3 3 0 0 1 1
*/
