import os
from subprocess import run

import numpy as np

import sciunit
from neuron import h
import neuron
from channelunit.capabilities import NModlChannel


F = 96485.33212  # C mol^-1
R = 8.314462618  # J mol^-1 K^-1


class ModelPatch(sciunit.Model, NModlChannel):
    nai = 10  # mM
    ki = 140  # mM
    cai = 100e-6  # mM

    def get_E_rev_name(self):
        ions = list(self.patch.psection()["ions"].keys())
        if self.ion_name not in ions:
            raise SystemError("Could not find %s in patch. I have only" % (ion_name)
                              + str(ions))
        return "e%s" % self.ion_name

    def get_E_rev_value(self, E_rev=None):
        if self.ion_name.lower() == "k":
            internal = self.ki
            valence = 1
        elif self.ion_name.lower() == "na":
            internal = self.nai
            valence = 1
        elif self.ion_name.lower() == "ca":
            internal = self.cai
            valence = 2
        external = self.external_conc
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
        conc_fact = np.log(self.external_conc/internal)
        E_rev = 1e3*R*(273.15+self.temperature)/(valence*F)*conc_fact
        return E_rev

    def compile_and_add(self, recompile):
        working_dir = os.getcwd()
        os.chdir(self.mod_path)
        if recompile:
            p = run('nrnivmodl')
        neuron.load_mechanisms(self.mod_path)
        os.chdir(working_dir)

    def __init__(self, path_to_mods, channel_name, ion_name,
                 external_conc=None, gbar_name="gbar", temp=22, recompile=True,
                 liquid_junction_pot=0, cvode=True, v_rest=-65, E_rev=None):
        """
        ion_name: str
            most common ions are: na (sodium), k, ca (sometimes Ca, if your model has different Ca pools). 
            if you specify nonspecific, you need to provide a value of E_rev
            It is important to pay attention whether the ion variable is specified with a lowercase
            or an uppercase letter, because the name of the reversal potential variable is constructed
            based on the ion name (na -> ena, k -> ek, ca -> eca, Ca -> eCa).
        """
        h.load_file("stdrun.hoc")
        self.dt = 0.01
        self.channel_name = channel_name
        self.mod_path = path_to_mods
        self.compile_and_add(recompile)
        self.patch = h.Section(name="patch")
        self.patch.L = 10
        self.patch.Ra = 100
        self.patch.diam = 10
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
        self.external_conc = external_conc
        if ion_name.lower() == "nonspecific":
            if isinstance(E_rev, int) or isinstance(E_rev, float):
                self.E_rev = E_rev
            else:
                raise SystemExit('Unable to proceed, if E_rev is unknown.')
        else:
            E_rev_name = self.get_E_rev_name()        
            self.E_rev = self.get_E_rev_value(E_rev)
            setattr(self.patch,  E_rev_name, self.E_rev)

        self.cvode = cvode
            

    def set_vclamp(self, dur1, v1, dur2, v2):
        self.vclamp.dur1 = dur1
        self.vclamp.amp1 = v1 - self.junction
        self.vclamp.dur2 = dur2
        self.vclamp.amp2 = v2 - self.junction

    def get_activation_steady_state(self, stimulation_levels: list,
                                    v_hold: float, t_stop:float,
                                    chord_conductance=False,
                                    duration=1000, sim_dt=0.001):
        """
        channel should be a neuron density mechanism, clamp a SEClamp neuron 
        object
        stimulation dict {amplitude: stim_duration}
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
        return self.normalize_to_one(max_current)

    def get_inactivation_steady_state(self, stimulation_levels: list,
                                      v_test: float, t_test:float,
                                      chord_conductance=False, sim_dt=0.001):
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
        return self.normalize_to_one(max_current)
                
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
