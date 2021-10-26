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
         
        return self.soma.psection()["ions"][ions[0]][name]

    def compile_and_add(self):
        working_dir = os.getcwd()
        os.chdir(mechanisms_path)
        p = run('nrnivmodl')
        neuron.load_mechanisms(mechanisms_path)
        os.chdir(working_dir)

    def __init__(self, path_to_mods, channel_name,
                 gbar_name="gbar", temp=22, E_rev=None):
        self.mod_path = path_to_mods
        self.compile_and_add()
        self.soma = h.Section(name="soma",
                              cell="my_patch")
        self.soma.L = 1
        self.soma.diam = 1
        self.soma.insert("extracellular")
        self.channel = self.soma.insert(channel_name)
        #set up channel conductance/permeability in case it is 0
        if self.soma.psection()["density_mechs"][channel_name][gbar] == 0:
            for seg in self.soma:
                from_mech = getattr(seg, channel_name)
                gbar_val = 0.001
                setattr(from_mech, gbar, gbar_val)
        self.temperature = temp
        self.vclamp = h.SEClamp(soma(0.5))
        
        E_rev_name = self.get_E_rev_name()        
        if E_rev is None:
            if E_rev_name is None:
                raise SystemExit('Unable to proceed, if E_rev is unknown.')
            else:
                val_E_rev = self.get_E_rev_val()
                if val_E_rev is not None:
                    self.E_rev = val_E_rev
                else:
                    raise SystemExit('Unable to proceed, if E_rev is unknown.')
        else:
            self.E_rev = E_rev
            if E_rev_name is not None:
                setattr(soma,  E_rev_name, E_rev)
            
                
