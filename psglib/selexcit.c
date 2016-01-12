/*   selexcit   -- selective excitation with presat option

	selective excitation is based on double PFG spin-echo
		(aka excitation sculpting)
			Shaka, et.al  JACS, 117, 4199 (1995)


Krish Krishnamurthy  (July 1995)  */


#include <standard.h>

pulsesequence()
{
	double		selpw,
			selpwr,
			hsgpwr,
			gzlvl1,
			gzlvl2,
			gstab,
			gt1,
			gt2;
	char	sspul[MAXSTR],
		satflg[MAXSTR],
		pwshape[MAXSTR];

	hsgpwr = getval("hsgpwr");
	getstr("sspul",sspul);
	selpw = getval("selpw");
	selpwr = getval("selpwr");
	gzlvl1 = getval("gzlvl1");
	gzlvl2 = getval("gzlvl2");
	gstab = getval("gstab");
	gt1 = getval("gt1");
	gt2 = getval("gt2");
	getstr("satflg",satflg);
	getstr("pwshape",pwshape);

	mod4(ct,oph);
	add(two,oph,v1);

status(A);

	delay(5.0e-5);
	if (sspul[0] == 'y')
	{
		zgradpulse(hsgpwr,1e-3);
		rgpulse(pw90,zero,rof1,rof1);
		zgradpulse(hsgpwr,1e-3);
	}

        delay(d1);
	if (satflg[0] == 'y')
	{
		obspower(satpwr);
		if (satfrq != tof)
		obsoffset(satfrq);
		rgpulse(satdly,zero,rof1,rof1);
		if (satfrq != tof)
		obsoffset(tof);
		obspower(tpwr);
	}


status(B);
	rgpulse(pw,oph,rof1,rof1);
        delay(gstab);
        zgradpulse(gzlvl1,gt1);
        delay(gstab-WFG_START_DELAY);
        obspower(selpwr);
        shaped_pulse(pwshape,selpw,oph,rof1,rof1);
        obspower(tpwr);
        delay(gstab-WFG_STOP_DELAY);
        zgradpulse(gzlvl1,gt1);
        delay(2*gstab);
        zgradpulse(gzlvl2,gt2);
        delay(gstab-WFG_START_DELAY);
        obspower(selpwr);
        shaped_pulse(pwshape,selpw,v1,rof1,rof2);
        obspower(tpwr);
        delay(gstab-WFG_STOP_DELAY);
        zgradpulse(gzlvl2,gt2);
        delay(gstab);
status(C);
}

