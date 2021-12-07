import os
from subprocess import run

import numpy as np

import sciunit
from neuron import h
import neuron
from channelunit.capabilities import NModlChannel


loc = os.path.dirname(os.path.abspath(__file__))
mechanisms_path = os.path.join(loc, 'mechanisms')


F = 96485.33212  # C mol^-1
R = 8.314462618  # J mol^-1 K^-1


class ModelPatch(sciunit.Model, NModlChannel):
    _nai = 10  # mM
    _ki = 140  # mM
    _cai = 100e-6  # mM

    def _find_E_rev_name(self):
        ions = list(self.patch.psection()["ions"].keys())
        if self.ion_name not in ions:
            raise SystemError("Could not find %s in patch. I have only" %
                              (ion_name) + str(ions))
        return "e%s" % self.ion_name

    def _find_E_rev_value(self, E_rev=None):
        if self.ion_name.lower() == "k":
            internal = self._ki
            valence = 1
        elif self.ion_name.lower() == "na":
            internal = self._nai
            valence = 1
        elif self.ion_name.lower() == "ca":
            internal = self._cai
            valence = 2
        elif E_rev is None:
            raise SystemExit("Unknown ion type %s. Only now na, k, ca and Ca"
                              % self.ion_name)
        elif self._external_conc is not None:
            raise SystemExit("Unknown ion type %s. Only now na, k, ca and Ca"
                             % self.ion_name)
        external = self._external_conc
        if external is None:
            if E_rev is None:
                name = "e%s" % self.ion_name
                return self.patch.psection()["ions"][self.ion_name][name][0]
            else:
                return E_rev
        elif not isinstance(external, int) and not isinstance(external, float):
            if E_rev is None:
                name = "e%s" % self.ion_name
                return self.patch.psection()["ions"][self.ion_name][name][0]
            else:
                return E_rev
        assert internal > 0
        # convert to mV
        conc_fact = np.log(self._external_conc/internal)
        E_rev = 1e3*R*(273.15+self.temperature)/(valence*F)*conc_fact
        return E_rev

    @property
    def ki(self):
        return self._ki

    @ki.setter
    def ki(self, value):
        self._ki = value
        if self.ion_name == "k" and self._external_conc is not None:
            self.E_rev = self._find_E_rev_value()

    @property
    def nai(self):
        return self._nai

    @nai.setter
    def nai(self, value):
        self._nai = value
        if self.ion_name == "na" and self._external_conc is not None:
            self.E_rev = self._find_E_rev_value()

    @property
    def cai(self):
        return self._cai

    @cai.setter
    def cai(self, value):
        self._cai = value
        if self.ion_name.lower() == "ca" and self._external_conc is not None:
            self.E_rev = self._find_E_rev_value()

    @property
    def Cai(self):
        return self._cai

    @cai.setter
    def Cai(self, value):
        self._cai = value
        if self.ion_name.lower() == "ca" and self._external_conc is not None:
            self.E_rev = self._find_E_rev_value()

    @property
    def external_conc(self):
        return self._external_conc

    @external_conc.setter
    def external_conc(self, value):
        self._external_conc = value
        self.E_rev = self._find_E_rev_value()
        
    def compile_and_add(self, path, recompile):
        working_dir = os.getcwd()
        os.chdir(path)
        if recompile:
            p = run('nrnivmodl')
        neuron.load_mechanisms(path)
        os.chdir(working_dir)

    def __init__(self, path_to_mods, channel_name, ion_name,
                 external_conc=None, gbar_name="gbar", temp=22, recompile=True,
                 liquid_junction_pot=0, cvode=True, v_rest=-65, E_rev=None):
        """
        ion_name: str
            most common ions are: na (sodium), k, ca (sometimes Ca, if 
            your model has different Ca pools). 
            if you specify nonspecific, you need to provide a value of E_rev
            It is important to pay attention whether the ion variable is 
            specified with a lowercase or an uppercase letter, because 
            the name of the reversal potential variable is constructed
            based on the ion name (na -> ena, k -> ek, ca -> eca, Ca -> eCa).
        """
        h.load_file("stdrun.hoc")
        self.dt = 0.01
        self.channel_name = channel_name
        self.mod_path = path_to_mods
        self.compile_and_add(self.mod_path, recompile)
        self.patch = h.Section(name="patch")
        self.patch.L = 1
        self.patch.Ra = 100
        self.patch.diam = 1
        self.patch.insert("pas")
        self.patch.e_pas = v_rest
        self.patch.g_pas = 1/60000
        self.channel = self.patch.insert(self.channel_name)
        self.junction = liquid_junction_pot
        self.base_directory = "validation_results"
        #set up channel conductance/permeability in case it is 0
        chan = self.patch.psection()["density_mechs"][channel_name]
        if gbar_name not in chan.keys():
            raise SystemExit('Unable to proceed, unknown %s conductance (gbar)'
                             % channel_name)

        if chan[gbar_name][0] == 0:
            for seg in self.patch:
                from_mech = getattr(seg, channel_name)
                gbar_val = 0.001
                setattr(from_mech, gbar_name, gbar_val)
        self.temperature = temp
        self.vclamp = h.SEClamp(self.patch(0.5))
        self.ion_name = ion_name
        self._external_conc = external_conc
        if ion_name.lower() == "nonspecific":
            if isinstance(E_rev, int) or isinstance(E_rev, float):
                self.E_rev = E_rev
            else:
                raise SystemExit('Unable to proceed, if E_rev is unknown.')
        else:
            E_rev_name = self._find_E_rev_name()
            self.E_rev = self._find_E_rev_value(E_rev)
            setattr(self.patch,  E_rev_name, self.E_rev)
        self.cvode = cvode
            

    def set_vclamp(self, dur1, v1, dur2, v2):
        self.vclamp.dur1 = dur1
        self.vclamp.amp1 = v1 - self.junction
        self.vclamp.dur2 = dur2
        self.vclamp.amp2 = v2 - self.junction

    def get_activation_steady_state(self, stimulation_levels: list,
                                    v_hold: float, t_stop:float,
                                    power: int,
                                    chord_conductance=False,
                                    duration=1000, sim_dt=0.001):
        """
        Function for running step experiments to determine steady-state
        activation curves.


        Diagram of the experiment:
              ________________________ stimulation_levels
             |________________________
             |________________________
             |________________________
             |________________________
             |________________________
        _____| <- v_hold

        simulation_levels: list
           list of voltages to test
        v_hold: float
          initial holding potential
        t_stop:
          stimulation duration
        power: int
          power coefficient of the activation gate.
          If experimental curves take into account the power of 
          the activation gate, then power is not 1
        chord_conductance: boolean
          in many experiments current is normalized by membrane voltage 
          minus the ion's reversal potential.
        duration: float
          duration of the simulation
        sim_dt: float
          for channels that can not be simulated using cvode. This value 
          should be small.
        """
        if self.cvode:
            h.cvode_active(1)
        else:
            h.cvode_active(0)
            h.dt = sim_dt
        h.celsius = self.temperature
        max_current = {}
        current = h.Vector()
        current.record(self.vclamp._ref_i, self.dt)
        time = h.Vector()
        time.record(h._ref_t, self.dt)
        delay = 200
        stim_start = int(delay/self.dt)
        for level in stimulation_levels:
            self.set_vclamp(delay, v_hold, duration, level)
            h.init()
            h.tstop = t_stop
            h.run()
            I = current.as_numpy()[stim_start:]
            out = self.extract_current(I, chord_conductance)
            max_current[level] = max(out)
        result = self.normalize_to_one(max_current)
        if power != 1:
            for key in result.keys():
                result[key] = result[key]**(1/power)
        return result

    def get_inactivation_steady_state(self, stimulation_levels: list,
                                      v_test: float, t_test:float,
                                      power: int,
                                      chord_conductance=False, sim_dt=0.001):
        """
        Function for running step experiments to determine steady-state
        inactivation curves.


        Diagram of the experiment
              ___ v_test
             |   |
             |   |
        _____|   |
        _____|   |
        _____|   |
        _____| <- stimulation_levels

        simulation_levels: list
           list of voltages to test
        v_test: float
          voltage level to test inactivated 
        t_test: float
          duration of test pulse
        power: int
          power coefficient of the inactivation gate.
          If experimental curves take into account the power 
          of the inactivation gate, then power is not 1
        chord_conductance: boolean
          in many experiments current is normalized by membrane voltage
          minus the ion's reversal potential.
        sim_dt: float
          for channels that can not be simulated using cvode. T
          his value should be small.
        """
        if self.cvode:
            h.cvode_active(1)
        else:
            h.cvode_active(0)
            h.dt = sim_dt

        h.celsius = self.temperature
        delay = 200
        stim_start = int(delay/self.dt)
        t_stop = delay +  t_test
        max_current = {}
        current = h.Vector()
        current.record(self.vclamp._ref_i, self.dt)
        time = h.Vector()
        time.record(h._ref_t, self.dt)
        for v_hold in stimulation_levels:
            self.set_vclamp(delay, v_hold, t_test, v_test)
            h.init()
            h.tstop = t_stop
            h.run()
            I = current.as_numpy()[stim_start:]
            out = self.extract_current(I, chord_conductance)
            max_current[v_hold] = max(out)
        result = self.normalize_to_one(max_current)
        if power != 1:
            for key in result.keys():
                result[key] = result[key]**(1/power)
        return result

                
    def extract_current(self, I, chord_conductance):
        if chord_conductance:
            return I/(self.vclamp.amp2 - self.E_rev)
        return abs(I)

    @staticmethod
    def normalize_to_one(current):
        factor = max(current.values())
        new_current = {}
        for key in current.keys():
            new_current[key] = current[key]/factor
        return new_current



class ModelPatchWithCa(ModelPatch):
    def __init__(self, path_to_mods, channel_name, ion_name,
                 external_conc=None, gbar_name="gbar", temp=22, recompile=True,
                 liquid_junction_pot=0, cvode=True, v_rest=-65):
        self.compile_and_add(mechanisms_path, True)
        super(ModelPatchWithCa, self).__init__(path_to_mods, channel_name,
                                                         ion_name, external_conc,
                                                         gbar_name, temp, recompile,
                                                         liquid_junction_pot, cvode,
                                                         v_rest, E_rev=None)
        if self.ion_name == "ca":
            self.patch.insert("cad")
            self.patch.cainf_cad = self._cai
            h.cao0_ca_ion = self._external_conc
        elif self.ion_name == "Ca":
            self.patch.insert("Cad")
            self.patch.Cainf_Cad = self._cai
            h.Cao0_Ca_ion = self._external_conc
        else:
            raise SystemExit("Unknown ion %s. I only know Ca and ca"
                             % self.ion_name )

    @property
    def cai(self):
        if self.ion_name == "ca":
            return self.patch.cainf_cad
        elif self.ion_name == "Ca":
            return self.patch.cainf_cad

    @cai.setter
    def cai(self, value):
        self._cai = value
        if self.ion_name.lower() == "ca":
            self.patch.cainf_cad = self._cai
        elif self.ion_name == "Ca":
            self.patch.cainf_Cad = self._cai

    @property
    def Cai(self):
        if self.ion_name == "ca":
            return self.patch.cainf_cad
        elif self.ion_name == "Ca":
            return self.patch.cainf_cad


    @cai.setter
    def Cai(self, value):
        self._cai = value
        if self.ion_name.lower() == "ca":
            self.patch.cainf_cad = self._cai
        elif self.ion_name == "Ca":
            self.patch.cainf_Cad = self._cai

    @property
    def external_conc(self):
        return  self._external_conc

    @external_conc.setter
    def external_conc(self, value):
        self._external_conc = value
        if self.ca_ion == "ca":
            h.cao0_ca_ion = self._external_conc
        elif self.ca_ion == "Ca":
            h.Cao0_Ca_ion = self._external_conc
        else:
            raise SystemExit("Unknown ion %s. I only know Ca and ca"
                             % self.ion_name )


class ModelCellAttachedPatch(ModelPatch):
    def __init__(self, path_to_mods, channel_name, ion_name,
                 external_conc=None, gbar_name="gbar", temp=22, recompile=True,
                 liquid_junction_pot=0, cvode=True, v_rest=-65, E_rev=None):
        super(ModelDendriteAttachedPatch, self).__init__(path_to_mods, channel_name,
                                                         ion_name, external_conc,
                                                         gbar_name, temp, recompile,
                                                         liquid_junction_pot, cvode,
                                                         v_rest, E_rev)



class ModelWholeCellPatch(ModelPatch):
    def __init__(self, path_to_mods, channel_name, ion_name,
                 external_conc=None, gbar_name="gbar", temp=22, recompile=True,
                 liquid_junction_pot=0, cvode=True, v_rest=-65, E_rev=None):

        super(ModelWholeCellPatch, self).__init__(path_to_mods, channel_name,
                                                  ion_name, external_conc, gbar_name,
                                                  temp, recompile, liquid_junction_pot,
                                                  cvode, v_rest, E_rev)
        self._L = 10
        self._diam = 10
        self.patch.L = self._L
        self.patch.Ra = 100
        self.patch.diam = self._diam

    @property
    def L(self):
        return self.patch.L

    @L.setter
    def L(self, value):
        self._L = value
        self.patch.L = self._L

    @property
    def diam(self):
        return self.patch.diam

    @diam.setter
    def diam(self, value):
        self._diam = value
        self.patch.diam = self._diam

    
