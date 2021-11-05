import os
from subprocess import run

import sciunit
from neuron import h
import neuron
from channelunit.capabilities import NModlChannel

class ModelPatch(sciunit.Model, NModlChannel):

    def get_E_rev_name(self):
        ions = list(self.soma.psection()["ions"].keys())
        if not len(ions):
            return None
        if len(ions) > 1:
            return None
        return "e%s" % ions[0]

    def get_E_rev_value(self):
        ions = list(self.soma.psection()["ions"].keys())
        try:
            name = "e%s" % ions[0]
        except IndexError:
            return None
        return self.soma.psection()["ions"][ions[0]][name][0]

    def compile_and_add(self, recompile):
        working_dir = os.getcwd()
        os.chdir(self.mod_path)
        if recompile:
            p = run('nrnivmodl')
        neuron.load_mechanisms(self.mod_path)
        os.chdir(working_dir)

    def __init__(self, path_to_mods, channel_name,
                 gbar_name="gbar", temp=22, E_rev=None, recompile=True,
                 liquid_junction_pot=10, cvode=True):
        """
        Liquid junction potential set to 10 unless otherwise specified
        """
        h.load_file("stdrun.hoc")
        self.dt = 0.01
        self.channel_name = channel_name
        self.mod_path = path_to_mods
        self.compile_and_add(recompile)
        self.soma = h.Section(name="soma")
        self.soma.L = 1
        self.soma.diam = 1
        self.channel = self.soma.insert(self.channel_name)
        self.junction = liquid_junction_pot
        self.base_directory = "validation_results"
        #set up channel conductance/permeability in case it is 0
        chan = self.soma.psection()["density_mechs"][channel_name]
        if gbar_name not in chan.keys():
            raise SystemExit('Unable to proceed, unknown %s conductance (gbar)'
                             % channel_name)

        if chan[gbar_name][0] == 0:
            for seg in self.soma:
                from_mech = getattr(seg, channel_name)
                gbar_val = 0.001
                setattr(from_mech, gbar_name, gbar_val)
        self.temperature = temp
        self.vclamp = h.SEClamp(self.soma(0.5))
        
        E_rev_name = self.get_E_rev_name()        
        if E_rev is None:
            if E_rev_name is None:
                raise SystemExit('Unable to proceed, if E_rev is unknown.')
            else:
                val_E_rev = self.get_E_rev_value()
                if val_E_rev is not None:
                    self.E_rev = val_E_rev
                else:
                    raise SystemExit('Unable to proceed, if E_rev is unknown.')
        else:
            self.E_rev = E_rev
            if E_rev_name is not None:
                setattr(self.soma,  E_rev_name, E_rev)

            self.cvode = cvode
            

    def set_vclamp(self, dur1, v1, dur2, v2):
        self.vclamp.dur1 = dur1
        self.vclamp.amp1 = v1 - self.junction
        self.vclamp.dur2 = dur2
        self.vclamp.amp2 = v2 - self.junction

    def get_activation_steady_state(self, stimulation_levels: list,
                                    v_init: float, t_stop:float,
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
            self.set_vclamp(delay, v_init, duration, level)
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
        for level in stimulation_levels:
            self.set_vclamp(delay, level, t_test, v_test)
            h.init()
            h.tstop = t_stop
            h.run()
            I = current.as_numpy()[stim_start:]
            out = self.extract_current(I, chord_conductance)
            max_current[level] = max(out)
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
