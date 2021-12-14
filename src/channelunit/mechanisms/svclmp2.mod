TITLE svclmp.mod
COMMENT
Single electrode Voltage clamp with TWELVE levels.
Clamp is on at time 0, and off at time
dur1+dur2+dur3+...+dur12. When clamp is off the injected current is 0.
The clamp levels are amp1, amp2, amp3, .., amp12.
i is the injected current, vc measures the control voltage)
Do not insert several instances of this model at the same location in order to
make level changes. That is equivalent to independent clamps and they will
have incompatible internal state values.
The electrical circuit for the clamp is exceedingly simple:
vc ---'\/\/`--- cell
        rs

Note that since this is an electrode current model v refers to the
internal potential which is equivalent to the membrane potential v when
there is no extracellular membrane mechanism present but is v+vext when
one is present.
Also since i is an electrode current,
positive values of i depolarize the cell. (Normally, positive membrane currents
are outward and thus hyperpolarize the cell)
ENDCOMMENT

INDEPENDENT {t FROM 0 TO 1 WITH 1 (ms)}

DEFINE NSTEP 3

NEURON {
	POINT_PROCESS SEClampOLS : with leak subtraction (Bezanilla and Clay 1977)
	ELECTRODE_CURRENT i
	RANGE dur1, amp1, dur2, amp2, dur3, amp3, rs, vc, i
	RANGE dur4, amp4, dur5, amp5, dur6, amp6, dur7, amp7
	RANGE dur8, amp8, dur9, amp9, dur10, amp10, dur11, amp11
	RANGE dur12, amp12
}

UNITS {
	(nA) = (nanoamp)
	(mV) = (millivolt)
	(uS) = (microsiemens)
}


PARAMETER {
	rs = 1 (megohm) <1e-9, 1e9>
	dur1 (ms) 	  amp1 (mV)
	dur2 (ms) <0,1e9> amp2 (mV)
	dur3 (ms) <0,1e9> amp3 (mV)
	dur4 (ms) <0,1e9> amp4 (mV)
	dur5 (ms) <0,1e9> amp5 (mV)
	dur6 (ms) <0,1e9> amp6 (mV)
	dur7 (ms) <0,1e9> amp7 (mV)
	dur8 (ms) <0,1e9> amp8 (mV)
	dur9 (ms) <0,1e9> amp9 (mV)
	dur10 (ms) <0,1e9> amp10 (mV)
	dur11 (ms) <0,1e9> amp11 (mV)
	dur12 (ms) <0,1e9> amp12 (mV)

}

ASSIGNED {
	v (mV)	: automatically v + vext when extracellular is present
	i (nA)
	vc (mV)
	tc2 (ms)
	tc3 (ms)
	tc4 (ms)
	tc5 (ms)
	tc6 (ms)
	tc7 (ms)
	tc8 (ms)
	tc9 (ms)
	tc10 (ms)
	tc11 (ms)
	tc12 (ms)
	on
}

INITIAL {
	tc2 = dur1 + dur2
	tc3 = tc2 + dur3
	tc4 = tc3 + dur4
	tc5 = tc4 + dur5
	tc6 = tc5 + dur6
	tc7 = tc6 + dur7
	tc8 = tc7 + dur8
	tc9 = tc8 + dur9
	tc10 = tc9 + dur10
	tc11 = tc10 + dur11
	tc12 = tc11 + dur12
	on = 0
}

BREAKPOINT {
	SOLVE icur METHOD after_cvode
	vstim()
}

PROCEDURE icur() {
	if (on) {
		i = (vc - v)/rs
	}else{
		i = 0
	}
}

COMMENT
The SOLVE of icur() in the BREAKPOINT block is necessary to compute
i=(vc - v(t))/rs instead of i=(vc - v(t-dt))/rs
This is important for time varying vc because the actual i used in
the implicit method is equivalent to (vc - v(t)/rs due to the
calculation of di/dv from the BREAKPOINT block.
The reason this works is because the SOLVE statement in the BREAKPOINT block
is executed after the membrane potential is advanced.

It is a shame that vstim has to be called twice but putting the call
in a SOLVE block would cause playing a Vector into vc to be off by one
time step.
ENDCOMMENT

PROCEDURE vstim() {
	on = 1
	if (dur1) {at_time(dur1)}
	if (dur2) {at_time(tc2)}
	if (dur3) {at_time(tc3)}
	if (dur4) {at_time(tc4)}
	if (dur5) {at_time(tc5)}
	if (dur6) {at_time(tc6)}
	if (dur7) {at_time(tc7)}
	if (dur8) {at_time(tc8)}
	if (dur9) {at_time(tc9)}
	if (dur10) {at_time(tc10)}
	if (dur11) {at_time(tc11)}
	if (dur12) {at_time(tc12)}
	if (t < dur1) {
		vc = amp1
	}else if (t < tc2) {
		vc = amp2
	}else if (t < tc3) {
		vc = amp3
	}else if (t < tc4) {
		vc = amp4
	}else if (t < tc5) {
		vc = amp5
	}else if (t < tc6) {
		vc = amp6
	}else if (t < tc7) {
		vc = amp7
	}else if (t < tc8) {
		vc = amp8
	}else if (t < tc9) {
		vc = amp9
	}else if (t < tc10) {
		vc = amp10
	}else if (t < tc11) {
		vc = amp11
	}else if (t < tc12) {
		vc = amp12
	}else {
		vc = 0
		on = 0
	}
	icur()
}

