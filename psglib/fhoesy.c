/* fhoesy - 2d fluorine proton NOE experiment
   corrected for double quantum artifacts
   Lix,B., Sonnichsen, F.D.,Sykes, B.D., J. Mag. Res.,
   121, series A, 1996, 83-87.  

   Parameters:

   pw:  19F pw90
   p1:  19F pw180
   d1:  relaxation delay
   d2:  evolution delay
   pwx: 1H pw90
   pwxlvl: x channel power level        
   mix:  NOESY mixing time
*/

#include <standard.h>

static int phs1[16] = {0, 3, 2, 1, 1, 0, 3, 2, 0, 3, 2, 1, 1, 0, 3, 2},
	   phs2[1]  = {1},
	   phs3[16] = {2, 3, 0, 1, 3, 0, 1, 2, 2, 3, 0, 1, 3, 0, 1, 2},
	   phs4[16] = {1, 1, 1, 1, 1, 1, 1, 1, 3, 3, 3, 3, 3, 3, 3, 3},
	   phs5[16] = {0, 0, 3, 3, 1, 1, 2, 2, 0, 0, 3, 3, 1, 1, 2, 2},
	   phs6[16] = {0, 2, 3, 1, 1, 3, 2, 0, 0, 2, 3, 1, 1, 3, 2, 0};

pulsesequence()

{
  double mix, pwx, pwxlvl, gzlvl1, gt1, gstab, corr;

  /* load variables */
  mix = getval("mix");
  pwx = getval("pwx");
  pwxlvl = getval("pwxlvl");
  gzlvl1 = getval("gzlvl1");
  gt1 = getval("gt1");
  gstab = getval("gstab");

  settable(t1, 16, phs1);
  settable(t2,  1, phs2);
  settable(t3, 16, phs3);
  settable(t4, 16, phs4);
  settable(t5, 16, phs5);
  settable(t6, 16, phs6);

  getelem(t1, ct, v1);
  getelem(t2, ct, v2);
  getelem(t3, ct, v3);
  getelem(t4, ct, v4);
  getelem(t5, ct, v5);
  getelem(t6, ct, oph);

  if (phase1 == 2)
    incr(v1);
  dbl(id2,v14);
  add(v1, v14, v1);
  add(oph, v14, oph);

  /* START PULSE SEQUENCE */
  status(A);
    delay(d1);
    decpower(pwxlvl);
  status(B);
    decrgpulse(pwx, v1, rof1, 0.0);
    corr = p1 / 2.0 + 2.0 * pwx / 3.14159 + rof1; 
    if ((d2 / 2.0 - corr) > 0.0)
      delay(d2 / 2.0 - corr);
    rgpulse(p1, v2, rof1, 0.0);
    if ((d2 / 2.0 - corr) > 0.0)
      delay(d2 / 2.0 - corr);
    decrgpulse(pwx, v3, rof1, 0.0);
    zgradpulse(gzlvl1, gt1);
  status(C);
    delay(mix);
    decrgpulse(pwx, v4, rof1, 0.0);
    rgpulse(pw, v5, rof1 + 0.0001, rof2);
    decpower(dpwr);
  status(D);
}
