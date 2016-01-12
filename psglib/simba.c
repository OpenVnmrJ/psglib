/* simba.c  - selincor & simba supported . F1 Region 2D Supported as well.
   TPPI implementations, steady-state transients for either all or
   just the first t1 increment, pulse+receiver phasecycling during
   steady-state transients, and presaturation using the transmitter.

   SELECTIVE OPTIONS via pwxs and refoc. SIMBA Supported RC 12 Sept91
   if pwxs = 0 then a normal hmqc/hmbc is done although the f1 phasing is not
   as good as with the hmqc.c sequence. WAVEFORM GENERATOR REQUIRED for a
   selective experiment. To date the best pulse is an EBURP although a
   UBURP would probably be better (Rolf K.). For selective HMBC (simba) try a
   refoc value of 100 (us) with a 6000 (us) pulse. For selective HMQC (selincor)
   set refoc = 3600 (us) if you are using an EBURP or 3600+1/2 the softpulse
   if you are using a gaussian.

   Trim X Trim Y sspulses preceed d1, recommended. Hypercomplex FAD is supported
   and highly recommended with either normal 2D or region-selected 2D.

   Ron Crouch (919) 248-3840.

   ref: Summers et al., j. amer. chem. soc. 108:4285-4294 (1986)
   J. Magn. Reson. 92:189-194(1991)

   Other than sstrim, pwxs, pwxslvl, and refoc, parameters are identical to
   hmqc.c

   1986-12-31 - s.l. patt
   1988-01-27 - revised (s.f.)
   1988-08-01 - revised
   1090-01-24 - revised
   1992-05    - revised for selincor and simba (R. Crouch)
   2003-07-01 - revised for VNMR 6.1 PSG (r.k.)
*/



#include <standard.h>

pulsesequence()
{
  /* VARIABLE DECLARATION */
  double at = getval("at"),
	 pwxlvl = getval("pwxlvl"),
	 pwx = getval("pwx"),
         pwxslvl = getval("pwxslvl"),
	 pwxs = getval("pwxs"),
	 refoc = getval("refoc"),
         j = getval("j"),
	 bird, total_time, decup_time,
	 null = getval("null"),
         satfrq = getval("satfrq"),
	 satdly = getval("satdly"),
	 satpwr = getval("satpwr"),
         sstrim = getval("sstrim"),
	 taumb = 0.0;
  int decup_ok;
  char satflg[MAXSTR], mbond[MAXSTR], shape[MAXSTR];

  getstr("satflg", satflg);
  getstr("mbond", mbond);
  getstr("shape", shape);

  /* INITIALIZE VARIABLES */
  if (mbond[0] == 'y')
    taumb = getval("taumb");

  if (j > 0.0)
  {
    bird = 1.0 / (2.0 * j);
  }
  else
  {
    bird = 0.0;
  }



  /* CHECK CONDITIONS */
  /* Check for correct system hardware configuration */
  if (!newdec)
  {
    printf("This sequence requires direct synthesis RF on DEC.\n");
    abort(1);
  }
  if (!newtransamp)
  {
    printf("This sequence requires linear amplifiers on XMTR.\n");
    abort(1);
  }

  /* Check for correct DM settings */
  if ((dm[A] == 'y') || (dm[B] == 'y'))
  {
    printf("DM must be set to either 'nny' or 'nnn'.\n");
    abort(1);
  }
  if ((mbond[0] == 'y') && (dm[C] == 'y'))
  {
    printf("DM must be set to 'nnn' for multiple-bond correlation.\n");
    abort(1);
  }

  /* Check for correct decoupler duty cycle */
  if (null != 0.0)
  {
    decup_time = 6.0 * pwx;
  }
  else
  {
    decup_time = 2.0 * pwx;
  }
  if (dm[C] == 'y')
    decup_time = decup_time + at;

  total_time = satdly + d1 + at + 2.0 * bird + null;
  if (null != 0.0)
    total_time = total_time + 2.0 * bird;

  decup_ok = (((decup_time / total_time) < 0.2) && (dpwr < 55.0));
  if (!decup_ok)
  {
    printf("Decoupler duty cycle must be under 20%%.\n");
    abort(1);
  }


  /* STEADY-STATE PHASECYCLING */
  sub(ct, ssctr, v10);
  hlv(v10, v1);
  dbl(v10, v3);


  /* PHASECYCLE CALCULATION */
  if (mbond[0] == 'y')
  {
    hlv(v1, v11);
    hlv(v11, v4);
    dbl(v11, v11);		/* J-filter phase subcycle */
  }
  else
  {
    hlv(v1, v4);
  }

  dbl(v4, v8);			/* suppression of 1H 180 artifacts */
  mod2(v1, v1);			/* QIS subcycle */

  hlv(v4, v5);
  hlv(v5, v5);
  dbl(v5, v5);			/* suppression of artifacts from 2nd X 90 */

  add(v1, one, v9);
  add(two, v1, v2);		/* for composite 180 */

  add(v5, v3, oph);
  add(v8, oph, oph);		/* receiver phasecycle */

  add(oph, v1, oph);
  add(v3, v1, v3);
  add(v4, v1, v4);
  add(v5, v1, v5);		/* addition of QIS subcycle */

  if (phase1 == 2)
    incr(v3);
  if (phase1 == 3)		/* TPPI with complex FT */
    add(v14, v3, v3);

  if ((phase1 == 1) || (phase1 == 2))	/*HYPER complex with FAD */
  {
    initval(2.0 * (double) ((int) (d2 * getval("sw1") + 0.5) % 2), v14);
    add(v3, v14, v3);
    add(oph, v14, oph);
  }


  /* BEGIN ACTUAL PULSE SEQUENCE CODE */
  status(A);
    decpower(pwxlvl);
    if (sstrim != 0.0)
    {
      obspower((tpwr - 4.0));
      rgpulse(sstrim, zero, 5.0e-5, 2e-6);
      rgpulse(sstrim, one, 2.0e-6, rof1);
    }
    obspower(tpwr);
    hsdelay(d1);

    /* selective saturation period */
    if (satflg[0] == 'y')
    {
      obsoffset(satfrq);
      obspower(satpwr);
      rgpulse(satdly, zero, 4.0e-5, 1.0e-5);
      obsoffset(tof);
      obspower(tpwr);
      delay(4.0e-5);
    }


  status(B);
    /* if null is 0 eliminate bird inversion pulse */
    /* bird pulse - 180 for protons not coupled to 13C */
    if ((null != 0.0) && (mbond[0] == 'n'))
    {
      rcvroff();
      rgpulse(pw, v1, rof1, 0.0);
      if (pw > pwx)
      {
        delay(bird - rof1 - 1.0e-6 - 1.5 * pw - pwx);
      }
      else
      {
        delay(bird - rof1 - 1.0e-6 - 2 * pwx - 0.5 * pw);
      }
      decrgpulse(pwx, v9, rof1, 0.0);
      simpulse(2.0 * pw, 2.0 * pwx, v1, v1, 1.0e-6, 0.0);
      decphase(v9);
      decrgpulse(pwx, v9, 1.0e-6, 0.0);
      if (pw > pwx)
      {
        delay(bird - rof1 - 1.0e-6 - 1.5 * pw - pwx);
      }
      else
      {
        delay(bird - rof1 - 1.0e-6 - 2 * pwx - 0.5 * pw);
      }
      rgpulse(pw, v2, rof1, rof2);
      rcvron();

      /* nulling time for protons not coupled to 13C */
      if (satflg[1] == 'y')
      {
        obsoffset(satfrq);
        obspower(satpwr);
        rgpulse(null, zero, 4.0e-5, 1.0e-5);
        obsoffset(tof);
        obspower(tpwr);
        delay(4.0e-5);
      }
      else
      {
        hsdelay(null);
      }
    }

    rcvroff();
    rgpulse(pw, v1, rof1, 0.0);
    delay(bird - 0.5 * pw - 0.5 * pwx - rof1);

    if (mbond[0] == 'y')		/* one-bond J(CH)-filter */
    {
      decrgpulse(pwx, v11, rof1, 0.0);
      delay(taumb - bird - rof1 - pwx);
    }

    txphase(v4);			/* presets transmitter phase */
    decrgpulse(pwx, v3, rof1, 0.0);
    if (pwxs != 0.0)
    {
      decpower(pwxslvl);		/* scale for low power */
    }
    delay(d2 / 2.0);
    if (pwxs != 0.0)
      delay(pwxs / 2.0);
    decphase(v5);
    rgpulse(2.0 * pw, v4, 0.0, 0.0);
    delay(d2 / 2.0);
    if (pwxs != 0.0)
    {
      decshaped_pulse(shape, pwxs, v5, 0.0, 0.0);
      delay(refoc);
    }
    else
    {
      decrgpulse(pwx, v5, 0.0, 0.0);
    }
    decpower(dpwr);
    rcvron();
    delay(rof1);
    if ((mbond[0] == 'n') && (pwxs == 0.0))
    {
      if (pwx > pw)
      {
        delay(bird - 4.2e-6 - 0.5 * pwx - rof2);
      }
      else
      {
        delay(bird - 4.2e-6 - 0.5 * pw - rof2);
      }
    }

  status(C);
}
