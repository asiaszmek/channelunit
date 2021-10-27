import os
from subprocess import run

import sciunit
from neuron import h
import neuron
import channelunit.capabilities as cap

class ModelPatch(sciunit.Model,
                 cap.NModlChannel):

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
        if chan[gbar_name] == 0:
            for seg in self.soma:
                from_mech = getattr(seg, channel_name)
                gbar_val = 0.001
                setattr(from_mech, gbar, gbar_val)
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
            
                
