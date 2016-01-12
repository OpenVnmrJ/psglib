/* fhoesy1d  fluorine proton 1D NOE experiment using steady-state NOE

   Parameters:

   pw:   19F pw90
   d1:   relaxation delay
   d2:   1H saturation pulse
   d3:   delay between saturation and 19F pulse (both pulses
         must go through high band amp) 
   dpwr: 1H channel saturation power level      
   dof:  frequency of 1H irradiation

*/




#include <standard.h>

static int phs1[4] = {0, 2, 1, 3},
           phs2[4] = {0, 2, 1, 3};



pulsesequence()
{
  double d2 = getval("d2"),
	 d3 = getval("d3");

  /* load variables */
  dof = getval("dof");
  dpwr = getval("dpwr");
  if (dpwr > 10)
  {
    printf("dpwr too large - acquisition aborted. \n");
    abort(1);
  }

  settable(t1, 4, phs1);
  settable(t2, 4, phs2);
  getelem(t1, ct, v1);
  getelem(t2, ct, oph);

  /* START PULSE SEQUENCE */
  status(A);
    delay(d1);
  status(B);
    decpower(dpwr);
    rcvroff();
    decon();
    delay(d2);
    decoff();
    delay(d3);
    rgpulse(pw, v1, rof1, rof2);
    rcvron();
  status(C);
}
