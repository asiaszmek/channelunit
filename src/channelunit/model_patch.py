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
                 gbar_name="gbar", temp=22, E_rev=None,
                 recompile=True):
        h.load_file("stdrun.hoc")
        self.mod_path = path_to_mods
        self.compile_and_add(recompile)
        self.soma = h.Section(name="soma")
        self.soma.L = 1
        self.soma.diam = 1
        self.soma.insert("extracellular")
        self.channel = self.soma.insert(channel_name)
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
        
    def get_activation_steady_state(self, stimulation_levels: list,
                                    v_init: float, t_stop:float,
                                    chord_conductance=False,
                                    duration=400, sim_dt=0.001):
        """
        channel should be a neuron density mechanism, clamp a SEClamp neuron 
        object
        stimulation dict {amplitude: stim_duration}
        """
        dt = 0.01
        h.celsius = self.temperature
        max_current = {}
        current = h.Vector().record(self.soma(0.5)._ref_i_membrane, dt)
        time = h.Vector().record(h._ref_t, dt)
        vm = h.Vector().record(self.soma(0.5)._ref_v, dt)
        delay = 200
        for level in stimulation_levels:
            self.vclamp.dur1 = delay
            self.vclamp.amp1 = v_init
            self.vclamp.dur2 = duration
            self.vclamp.amp2 = level
            h.dt = sim_dt
            h.tstop = t_stop
            h.run(t_stop)
            if chord_conductance:
                out = current.as_numpy()/(vm - self.E_rev)
            else:
                out = abs(current.as_numpy())
            max_current[level] = max(out)
        return max_current

    def get_inactivation_steady_state(self, stimulation_levels: list,
                                      v_test: float, t_test:float,
                                      chord_conductance=False, sim_dt=0.001):
        h.celsius = self.temperature
        dt = 0.01
        delay = 200
        t_stop = delay +  t_test
        max_current = {}
        current = h.Vector().record(self.soma(0.5)._ref_i_membrane, dt)
        time = h.Vector().record(h._ref_t, dt)
        vm = h.Vector().record(self.soma(0.5)._ref_v, dt)
        for level in stimulation_levels:
            self.vclamp.dur1 = delay
            self.vclamp.amp1 = level
            self.vclamp.dur2 = t_test
            self.vclamp.amp2 = v_test
            h.dt = sim_dt
            h.tstop = t_stop
            h.run(t_stop)
            if chord_conductance:
                out = current.as_numpy()/(vm - self.E_rev)
            else:
                out = abs(current.as_numpy())
            max_current[level] = max(out)
        return max_current
                
