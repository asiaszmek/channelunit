TITLE Calcium decay
: as described in Bhalla and Bower, J. Neurophysiol. 69:1948-1983 (1993)
: written by Andrew Davison
: partially based on cadecay.mod by Alain Destexhe, Salk Institute 1995.
: 25-08-98

INDEPENDENT {t FROM 0 TO 1 WITH 1 (ms) }

NEURON{
	SUFFIX Cad
	USEION Ca READ iCa, Cai WRITE Cai VALENCE 2
	RANGE iCa, channel_flow, depth, B
	RANGE Cai, tau, Cainf
}

UNITS {
	(mA) = (milliamp)
	(mV) = (millivolt)
	(molar) = (1/liter)
	(mM) = (millimolar)
	(um) = (micron)
}

CONSTANT {
        FARADAY = 96154 (coul)
	:FARADAY = 93149 (coul)		: moles do not appear in units
					: note this value is chosen to fit with
					: Genesis
}

PARAMETER {
	dt (ms)
	depth = 1 	(um)		: shell within which Cai is Calculated
					: to match Bhalla and Bower 1993 set
					: depth = diam/4 for each compartment
	tau = 100 	(ms)		: Cai deCay constant
	Cainf = 100e-6	(mM)		: baseline Calcium concentration
	iCa		(mA/cm2)
}

STATE {
	Cai		(mM)
}

INITIAL {
	Cai = Cainf
}

ASSIGNED {
	channel_flow	(mM/ms)
	B		(mM cm2/ms/mA)
}

BREAKPOINT {
	SOLVE state METHOD cnexp
}

DERIVATIVE state {
	B = -(1e4)/(2*FARADAY*depth)
	channel_flow = B*iCa
	if (channel_flow <= 0.0 ) { channel_flow = 0.0 }	: one way flow in channel
	Cai' = channel_flow  - (Cai - Cainf)/tau
}
	











