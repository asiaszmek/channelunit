TITLE Low threshold calcium current
: Cav3.1?
:   Ca++ current responsible for low threshold spikes (LTS)
:   Differential equations
:
:   Model of Huguenard & McCormick, J Neurophysiol 68: 1373-1383, 1992.
:   The kinetics is described by Goldman-Hodgkin-Katz equations,
:   using a m2h format, according to the voltage-clamp data
:   (whole cell patch clamp) of Huguenard & Prince, J. Neurosci. 
:   12: 3804-3817, 1992.
:
:   This model is described in detail in:
:   Destexhe A, Neubig M, Ulrich D and Huguenard JR.  
:   Dendritic low-threshold calcium currents in thalamic relay cells.  
:   Journal of Neuroscience 18: 3574-3588, 1998.
:   (a postscript version of this paper, including figures, is available on
:   the Internet at http://cns.fmed.ulaval.ca)
:
:    - shift parameter for screening charge
:    - empirical correction for contamination by inactivation (Huguenard)
:    - GHK equations
:
:
:   Written by Alain Destexhe, Laval University, 1995
:

: 2019: From ModelDB, accession no. 279
: Modified qm and qh by Elisabetta Iavarone @ Blue Brain Project
: See PARAMETER section for references

UNITS {
	(molar) = (1/liter)
	(mV) =	(millivolt)
	(mA) =	(milliamp)
	(mM) =	(millimolar)

	FARADAY = (faraday) (coulomb)
	R = (k-mole) (joule/degC)
}


NEURON {
        THREADSAFE
	SUFFIX caT
	USEION ca READ cai, cao WRITE ica VALENCE 2
	RANGE gbar, ica

}


PARAMETER {
	gbar	=.2e-3	(cm/s)	: Maximum Permeability
	shift	= 2 	(mV)	: corresponds to 2mM ext Ca++
	actshift = 0 	(mV)	: shift of activation curve (towards hyperpol)
	qm      = 2.5		: Amarillo et al., J Neurophysiol, 2014
	qh      = 2.5           : Amarillo et al., J Neurophysiol, 2014
        km      = 6.2
	kh      = 4.0
        T_denom = 10 (degC)
        T_thresh = 23.5 (degC)					  
}

ASSIGNED {
        v (mV)
	ica	(mA/cm2)
	celsius (degC)
        cai (mM)
        cao (mM)
	m_inf
	tau_m	(ms)
	h_inf
	tau_h	(ms)
	phi_m
	phi_h
}



STATE {
	m h
}



BREAKPOINT {
	SOLVE castate METHOD cnexp
	ica = gbar * m*m*h * ghk(v, cai, cao)
}

INITIAL {
     	phi_m = qm ^ ((celsius-T_thresh)/T_denom)
	phi_h = qh ^ ((celsius-T_thresh)/T_denom)
        evaluate_fct()
	m = m_inf
	h = h_inf
}

DERIVATIVE castate {
	evaluate_fct()
	m' = (m_inf - m) / tau_m
	h' = (h_inf - h) / tau_h
}




UNITSOFF

PROCEDURE evaluate_fct() {
:
:   The kinetic functions are taken as described in the model of 
:   Huguenard & McCormick, and corresponds to a temperature of 23-25 deg.
:   Transformation to 36 deg assuming Q10 of 5 and 3 for m and h
:   (as in Coulter et al., J Physiol 414: 587, 1989).
:
:   The activation functions were estimated by John Huguenard.
:   The V_1/2 were of -57 and -81 in the vclamp simulations, 
:   and -60 and -84 in the current clamp simulations.
:
:   The activation function were empirically corrected in order to account
:   for the contamination of inactivation.  Therefore the simulations 
:   using these values reproduce more closely the voltage clamp experiments.
:   (cfr. Huguenard & McCormick, J Neurophysiol, 1992).
:
	m_inf = 1.0 / ( 1 + exp(-(v+shift+actshift+57)/km) )
	h_inf = 1.0 / ( 1 + exp((v+shift+81)/kh) )

	tau_m = ( 0.612 + 1.0 / ( exp(-(v+shift+actshift+132 (mV))/16.7 (mV)) + exp((v+shift+actshift+16.8 (mV))/18.2 (mV)) ) ) / phi_m
	if( (v+shift) < -80 (mV)) {
		tau_h = exp((v+shift+467 (mV))/66.6 (mV)) / phi_h
	} else {
		tau_h = ( 28 (mV) + exp(-(v+shift+22 (mV))/10.5 (mV)) ) / phi_h
	}

	: EI compare with tau_h on ModelDB, no. 3817
}
UNITSON
		     
FUNCTION ghk(v(mV), ci(mM), co(mM)) (.001 coul/cm3) {
	LOCAL z, eci, eco
	z = (1e-3)*2*FARADAY*v/(R*(celsius+273.15))
	eco = co*efun(z)
	eci = ci*efun(-z)
	:high cao charge moves inward
	:negative potential charge moves inward
	ghk = (.001)*2*FARADAY*(eci - eco)
}

FUNCTION efun(z) {
	if (fabs(z) < 1e-4) {
		efun = 1 - z/2
	}else{
		efun = z/(exp(z) - 1)
	}
}

