TITLE HVA L-type calcium current (Cav1.2)


UNITS {
    (mV) = (millivolt)
    (mA) = (milliamp)
    (S) = (siemens)
    (molar) = (1/liter)
    (mM) = (millimolar)
    FARADAY = (faraday) (coulomb)
    R = (k-mole) (joule/degC)
}

NEURON {
    THREADSAFE
    SUFFIX TCcal
    USEION ca READ cai, cao WRITE ica VALENCE 2
    RANGE gbar, ica
    
}

PARAMETER {
    gbar = 1e-4 (cm/s)
    T_denom = 10 (degC)
} 

ASSIGNED { 
    v (mV)
    ica (mA/cm2)
    celsius (degC)
    cai (mM)
    cao (mM)
    minf
    mtau (ms)
    tcorr

}

STATE { m  }

BREAKPOINT {
    SOLVE states METHOD cnexp
    ica = gbar*m*m*ghk(v, cai, cao)
}

INITIAL {
    rates()
    m = minf
    
}

DERIVATIVE states { 
    rates()
    m' = (minf-m)/mtau
   
}


UNITSOFF
PROCEDURE rates() {LOCAL a, b
    tcorr = 3^((celsius-23.5)/T_denom)			 
    a = 1.6 / (1+ exp(-0.072*(v-5)))
    b = 0.02 * vtrap( -(v-1.31), 5.36)

    mtau = 1/(a+b) / tcorr
    minf = 1/(1+exp((v+10)/-10))
}


FUNCTION vtrap(x,c) { 
	: Traps for 0 in denominator of rate equations
        if (fabs(x/c) < 1e-6) {
          vtrap = c + x/2 }
        else {
          vtrap = x / (1-exp(-x/c)) }
}
 

UNITSON	

FUNCTION ghk(v (mV), ci (mM), co (mM)) (.001 coul/cm3) {
    LOCAL z, eci, eco
    z = (1e-3)*2*FARADAY*v/(R*(celsius+273.15))
    if(z == 0) {
        z = z+1e-6
    }
    eco = co*(z)/(exp(z)-1)
    eci = ci*(-z)/(exp(-z)-1)
    ghk = (1e-3)*2*FARADAY*(eci-eco)
}



