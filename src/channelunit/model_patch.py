import os
from subprocess import run

import numpy as np

import sciunit
from neuron import h
from neuron import rxd
import neuron
from channelunit.capabilities import NModlChannel

N = 4
DT = 0.1
membrane_shell_width = .1

loc = os.path.dirname(os.path.abspath(__file__))
mechanisms_path = os.path.join(loc, 'mechanisms')

F = 96485.33212  # C mol^-1
R = 8.314462618  # J mol^-1 K^-1

class ModelPatch(sciunit.Model):
    def __init__(self, temp=22, R_m=20000,
                 v_rest=-65, liquid_junction_pot=0,
                 cvode=True,
                 sim_dt=0.001):
        h.load_file("stdrun.hoc")
        self.dt = DT
        self.compile_and_add(mechanisms_path, True)
        self.junction = liquid_junction_pot
        self.patch = h.Section(name="patch")
        self.patch.L = 1
        self.patch.Ra = 100
        self.patch.diam = 1
        self.patch.insert("pas")
        self.patch.e_pas = v_rest
        self.patch.g_pas = 1/R_m
        self.temperature = temp
        self.v_low = 29.5
        self.vclamp = h.SEClampOLS(self.patch(0.5))
        self.cvode = cvode
        if cvode:
            h.CVode()
            h.CVode().atol(1e-7)
        else:
            h.dt = sim_dt
        
    def compile_and_add(self, path, recompile):
        working_dir = os.getcwd()
        os.chdir(path)
        if recompile:
            p = run('nrnivmodl')
        neuron.load_mechanisms(path)
        os.chdir(working_dir)

    def set_vclamp(self, dur1, v1, dur2, v2, leak_subtraction, delay=200):
        self.vclamp.dur1 = dur1
        self.vclamp.amp1 = v1 - self.junction
        self.vclamp.dur2 = dur2
        self.vclamp.amp2 = v2 - self.junction
        if not leak_subtraction:
            self.vclamp.dur3 = 0
            self.vclamp.dur4 = 0
            self.vclamp.dur5 = 0
            self.vclamp.dur6 = 0
            self.vclamp.dur7 = 0
            self.vclamp.dur8 = 0
            self.vclamp.dur9 = 0
            self.vclamp.dur10 = 0
            self.vclamp.dur11 = 0
            self.vclamp.dur12 = 0
            return dur1+dur2+delay

        pulse_amp = self.vclamp.amp2 - self.vclamp.amp1
        self.vclamp.dur3 = delay
        self.vclamp.amp3 = self.vclamp.amp1
        v_sub = self.vclamp.amp1 - pulse_amp/4 - self.v_low 
        v_pulse = v_sub + pulse_amp/4 
        self.vclamp.amp4 = v_sub
        self.vclamp.dur4 = delay
        # 1st pulse
        self.vclamp.amp5 = v_pulse
        self.vclamp.dur5 = dur2
        self.vclamp.amp6 = v_sub
        self.vclamp.dur6 = delay
        #2nd pulse
        self.vclamp.amp7 = v_pulse
        self.vclamp.dur7 = dur2
        self.vclamp.amp8 = v_sub
        self.vclamp.dur8 = delay
        #3rd pulse
        self.vclamp.amp9 = v_pulse
        self.vclamp.dur9 = dur2
        self.vclamp.amp10 = v_sub
        self.vclamp.dur10 = delay
        #4th pulse
        self.vclamp.amp11 = v_pulse
        self.vclamp.dur11 = dur2
        self.vclamp.amp12 = v_sub
        self.vclamp.dur12 = delay
        return dur1 + 6*delay + 5*dur2

    def run(self, t_stop, dt=DT):
        current = h.Vector()
        current.record(self.vclamp._ref_i, dt)
        if self.cvode:
            h.CVode().re_init()
        else:
            h.init()
        h.tstop = t_stop
        h.run()
        return current.as_numpy()

    @classmethod
    def curr_stim_response(self, I, dur1, dur2, dt):
        current = I[int(dur1/dt)+10:int((dur1+dur2)/dt)]
        return current.copy()

    @classmethod
    def curr_leak_amp(self, I, dur1, dur2, delay, dt):
        t_start = dur1 + dur2 + 2*delay
        length = int((t_start+dur2)/dt) - int(t_start/dt) -10
        pulse = np.zeros((length,))
        for i in range(N):
            pulse += I[int(t_start/dt)+10:int((t_start+dur2)/dt)].copy()
            t_start += dur2 + delay
        return pulse


class ModelPatchWithChannel(ModelPatch, NModlChannel):
    _nai = 10  # mM
    _ki = 140  # mM
    _cai = 100e-6  # mM
    def __init__(self, path_to_mods: str, channel_name: str, ion_name: str,
                 external_conc=None, gbar_name="gbar", temp=22,
                 recompile=True, liquid_junction_pot=0,
                 cvode=True, R_m=20000, v_rest=-65, gbar_value=None,
                 directory="validation_results",  sim_dt=0.001):
        """
        ion_name: str
            most common ions are: na (sodium), k, ca (sometimes Ca, if 
            your model has different Ca pools). 
            if you specify nonspecific, you need to provide a value of E_rev
            It is important to pay attention whether the ion variable is 
            specified with a lowercase or an uppercase letter, because 
            the name of the reversal potential variable is constructed
            based on the ion name (na -> ena, k -> ek, ca -> eca, Ca -> eCa).
        R_m: float 
            mebrane resistivity (in ohm*cm^2)
        """

        self.channel_name = [channel_name]
        self.mod_path = path_to_mods
        self.compile_and_add(self.mod_path, recompile)
        super(ModelPatchWithChannel, self).__init__(temp=temp, 
                                                    R_m=R_m, v_rest=v_rest,
                                                    liquid_junction_pot=liquid_junction_pot,
                                                    cvode=cvode,
                                                    sim_dt=sim_dt,)

        self.channel = self.patch.insert(channel_name)
        self.base_directory = directory
        #set up channel conductance/permeability in case it is 0
        chan = self.patch.psection()["density_mechs"][channel_name]
        self.gbar_name = {channel_name: gbar_name}
        self.ion_name = ion_name
        self._external_conc = external_conc
        if gbar_name not in chan.keys():
            raise SystemExit('Unable to proceed, unknown %s conductance (gbar)'
                             % channel_name)

        if chan[gbar_name][0] == 0 or gbar_value is not None:
            for seg in self.patch:
                from_mech = getattr(seg, channel_name)
                if gbar_value is None:
                    gbar_value = 0.001
                setattr(from_mech, gbar_name, gbar_value)
        self.ca = None
        self.E_Ca_rev = None
        self._Ca_ext = None
        
    @property
    def gbar(self, channel_name=None):
        if channel_name is None:
            channel_name=self.channel_name[0]
        mech = self.patch.psection()["density_mechs"][channel_name]
        values = mech[self.gbar_name[channel_name]]
        if len(values) == 1:
            return values[0]
        elif len(set(values)) == 1:
            return values[0]
        return values

    @gbar.setter
    def gbar(self, value, channel_name=None):
        if channel_name is None:
            channel_name=self.channel_name[0]
        if not isinstance(value, list):
            value = [value]
        if len(value) == self.patch.nseg:
            for i, seg in enumerate(self.patch):
                from_mech = getattr(seg, channel_name)
                setattr(from_mech, self.gbar_name[channel_name], value[i])
        elif len(value) == 1:
            for i, seg in enumerate(self.patch):
                from_mech = getattr(seg, channel_name)
                setattr(from_mech, self.gbar_name[channel_name], value[0])
       
    def _find_E_rev_name(self):
        ions = list(self.patch.psection()["ions"].keys())
        if self.ion_name not in ions:
            raise SystemError("Could not find %s in patch. I have only" %
                              (ion_name) + str(ions))
        return "e%s" % self.ion_name

    def calc_E_rev(self, E_rev=None, external=None):
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
        
        if external is None:
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
        conc_fact = np.log(external/internal)
        E_rev = 1e3*R*(273.15+self.temperature)/(valence*F)*conc_fact
        return E_rev

    def generate_fname(self, suffix, stim_beg, stim_end, chord_conductance,
                       electrode_current, ca_conc):
        fname = suffix
        if chord_conductance:
            fname = "%s_%s" % (fname, "g")
        else:
            fname = "%s_%s" % (fname, "I")
        if not electrode_current:
            fname = "%s_%s" % (fname, "np_o")
        for channel_name in self.channel_name:
            fname = "%s_%s" % (fname, channel_name)
        if self._external_conc is not None:
            fname = "%s_%s_ext_%4.2f_mM" % (fname, self.ion_name,
                                        self._external_conc)
        elif self.E_rev is not None:
            fname = "%s_E_%s_rev_%4.2f_mV" % (fname, self.ion_name,
                                              self.E_rev)
        if self._Ca_ext is not None:
            fname = "%s_Ca_ext_%4.2f_mM" % (fname, self._Ca_ext)
        elif self.E_Ca_rev is not None:
            fname = "%s_E_Ca_rev_%4.2f_mV" % (fname, self._Ca_ext)

        
        fname = "%s_from_%4.2f_mV_to_%4.2f_mV" % (fname, stim_beg, stim_end)
        
        if self.cvode:
            fname = "%s_cvode" % fname
        if ca_conc:
            fname = "%s_Ca_conc" % fname
        return fname

    def get_activation_traces(self, stimulation_levels: list,
                              v_hold: float, t_stop:float,
                              chord_conductance=False,
                              electrode_current=True,
                              interval=200,
                              save_traces=True, save_ca=True):
        """
        Function for running step experiments to determine 
        current/chord conductance traces


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
        chord_conductance: boolean
          in many experiments current is normalized by membrane voltage 
          minus the ion's reversal potential.
        duration: float
          duration of the simulation
        """
        if save_traces:
            fname = self.generate_fname("Activation_traces",
                                        min(stimulation_levels),
                                        max(stimulation_levels),
                                        chord_conductance,
                                        electrode_current, False)
            output = []
            header = "time"
        if save_ca and self.ca is None:
            save_ca = False
        if save_ca:
            ca_fname = self.generate_fname("Activation_traces",
                                        min(stimulation_levels),
                                        max(stimulation_levels),
                                        chord_conductance,
                                        electrode_current, True)
            calcium_vals = []
            calcium = h.Vector()
            calcium.record(self.ca[self.memb_shell].nodes[0]._ref_concentration,
                           self.dt)
        else:
            ca_fname = ""
        h.celsius = self.temperature
        current_vals = {}
        current = h.Vector()
        if electrode_current:
            current.record(self.vclamp._ref_i, self.dt)
            leak_subtraction = True
        else:
            leak_subtraction = False
            if self.ion_name == "na":
                current.record(self.patch(0.5)._ref_ina, self.dt)
            elif self.ion_name == "k":
                current.record(self.patch(0.5)._ref_ik, self.dt)
            elif self.ion_name == "ca":
                current.record(self.patch(0.5)._ref_ica, self.dt)
            elif self.ion_name == "Ca":
                current.record(self.patch(0.5)._ref_iCa, self.dt)
            elif self.ion_name == "ba":
                current.record(self.patch(0.5)._ref_ica, self.dt)
            elif self.ion_name == "Ba":
                current.record(self.patch(0.5)._ref_iCa, self.dt)
            else:
                current.record(self.vclamp._ref_i, self.dt)
                leak_subtraction = True
        time = h.Vector()
        time.record(h._ref_t, self.dt)
        delay = 200
        stim_start = int(delay/self.dt)
        for i, level in enumerate(stimulation_levels):
            stim_stop = self.set_vclamp(delay, v_hold, t_stop, level,
                                        leak_subtraction,
                                        delay=interval)
            if self.cvode:
                h.CVode().re_init()
            else:
                h.init()
            h.tstop = stim_stop
            h.run()
            I = current.as_numpy()
            out = self.extract_current(I, chord_conductance, leak_subtraction,
                                       delay, t_stop, interval, self.dt)
            beg = int(delay/self.dt)+10
            end = int((delay+t_stop)/self.dt)
            if save_ca:
                if not i:
                    save_time = time.as_numpy()[beg: end].copy()
                    calcium_vals.append(save_time)
  
                calcium_vals.append(calcium.as_numpy()[beg: end].copy())
            if save_traces:
                if not i:
                    save_time = time.as_numpy()[beg:end].copy()                  
                    output.append(save_time)
                output.append(out)
                header += ";%4.2f" % level
            current_vals[level] = out
        if save_traces:
            path = os.path.join(self.base_directory, "data")
            if not os.path.exists(path):
                os.makedirs(path)
            path_to_save = os.path.join(path, "%s.csv" % fname)
            np.savetxt(path_to_save, np.array(output), delimiter=";",
                       header=header, comments="")
        if save_ca:
            path = os.path.join(self.base_directory, "data")
            if not os.path.exists(path):
                os.makedirs(path)
            path_to_save = os.path.join(path, "%s.csv" % ca_fname)
            np.savetxt(path_to_save, np.array(calcium_vals), delimiter=";",
                       header=header, comments="")
            
        return current_vals

    def get_activation_steady_state(self, stimulation_levels: list,
                                    v_hold: float, t_stop:float,
                                    power: int, chord_conductance,
                                    electrode_current, interval=200,
                                    normalization="to_one",
                                    save_traces=True, save_ca=True):
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
          should be small.
        """
        currents = self.get_activation_traces(stimulation_levels,
                                              v_hold, t_stop,
                                              chord_conductance,
                                              electrode_current,
                                              interval=interval,
                                              save_traces=save_traces,
                                              save_ca=save_ca)
        max_current = self.get_max_of_dict(currents, self.ion_name,
                                           chord_conductance)
        result = self.normalize_to_one(max_current, normalization)
        if power != 1:
            for key in result.keys():
                result[key] = result[key]**(1/power)
        return result

    def get_inactivation_traces(self, stimulation_levels: list,
                                v_test: float, t_test:float,
                                chord_conductance,
                                electrode_current,
                                interval=200,
                                save_traces=True, save_ca=True):
        """
        Function for running step experiments to determine steady-state
        inactivation currents.


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
        chord_conductance: boolean
          in many experiments current is normalized by membrane voltage
          minus the ion's reversal potential.
        """
        if save_traces:
            fname = self.generate_fname("Inactivation_traces",
                                        min(stimulation_levels),
                                        max(stimulation_levels),
                                        chord_conductance,
                                        electrode_current, False)
            output = []
            header = "time"
        if save_ca and self.ca is None:
            save_ca = False
        if save_ca:
            ca_fname = self.generate_fname("Inactivation_traces",
                                           min(stimulation_levels),
                                           max(stimulation_levels),
                                           chord_conductance,
                                           electrode_current, True)
            calcium_vals = []
            calcium = h.Vector()
            calcium.record(self.ca[self.memb_shell].nodes[0]._ref_concentration,
                           self.dt)
        else:
            ca_fname = ""

        h.celsius = self.temperature
        delay = 200
        stim_start = int(delay/self.dt)
        current_values = {}
        current = h.Vector()
        if electrode_current:
            current.record(self.vclamp._ref_i, self.dt)
            leak_subtraction = True
        else:
            leak_subtraction = False
            if self.ion_name == "na":
                current.record(self.patch(0.5)._ref_ina, self.dt)
            elif self.ion_name == "k":
                current.record(self.patch(0.5)._ref_ik, self.dt)
            elif self.ion_name == "ca":
                current.record(self.patch(0.5)._ref_ica, self.dt)
            elif self.ion_name == "Ca":
                current.record(self.patch(0.5)._ref_iCa, self.dt)
            elif self.ion_name == "ba":
                current.record(self.patch(0.5)._ref_ica, self.dt)
            elif self.ion_name == "Ba":
                current.record(self.patch(0.5)._ref_iCa, self.dt)
            else:
                current.record(self.vclamp._ref_i, self.dt)
                leak_subtraction = True
        time = h.Vector()
        time.record(h._ref_t, self.dt)
        for i, v_hold in enumerate(stimulation_levels):
            t_stop = self.set_vclamp(delay, v_hold, t_test, v_test,
                                     leak_subtraction,
                                     delay=interval)
            if self.cvode:
                h.CVode().re_init()
            else:
                h.init()
            h.tstop = t_stop
            h.run()
            I = current.as_numpy()
            out = self.extract_current(I, chord_conductance,
                                       leak_subtraction, delay, t_test,
                                       0, self.dt)
            current_values[v_hold] = out
            beg = int(delay/self.dt)+10
            end = int((delay+t_test)/self.dt)
            if save_ca:
                if not i:
                    save_time = time.as_numpy()[beg: end].copy()
                    calcium_vals.append(save_time)
                calcium_vals.append(calcium.as_numpy()[beg: end].copy())
            if save_traces:
                if not i:
                    save_time = time.as_numpy()[beg:end].copy()
                    output.append(save_time)
                output.append(out)
                header += ";%4.2f" % v_hold

        if save_traces:
            path = os.path.join(self.base_directory, "data")
            if not os.path.exists(path):
                os.makedirs(path)
            path_to_save = os.path.join(path, "%s.csv" % fname)
            np.savetxt(path_to_save, np.array(output), delimiter=";",
                       header=header, comments="")
        if save_ca:
            path = os.path.join(self.base_directory, "data")
            if not os.path.exists(path):
                os.makedirs(path)
            path_to_save = os.path.join(path, "%s.csv" % ca_fname)
            np.savetxt(path_to_save, np.array(calcium_vals), delimiter=";",
                       header=header, comments="")
  
        return current_values
        
    def get_inactivation_steady_state(self, stimulation_levels: list,
                                      v_test: float, t_test:float,
                                      power: int,
                                      chord_conductance=False,
                                      leak_subtraction=True,
                                      electrode_current=True,
                                      interval=200,
                                      normalization="to_one",
                                      save_traces=True, save_ca=True):
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
          his value should be small.
        """
        currents = self.get_inactivation_traces(stimulation_levels,
                                                v_test, t_test,
                                                chord_conductance,
                                                leak_subtraction,
                                                interval,
                                                save_traces=save_traces,
                                                save_ca=save_ca)
        max_current = self.get_max_of_dict(currents, self.ion_name,
                                           chord_conductance)
        result = self.normalize_to_one(max_current, normalization)
        if power != 1:
            for key in result.keys():
                result[key] = result[key]**(1/power)
        return result
                
    def extract_current(self, I, chord_conductance, leak_subtraction, dur1,
                        dur2, delay, dt):
        #dt = self.dt
        current = self.curr_stim_response(I, dur1, dur2, dt)
        
        #either step injection or the short pulse
        if leak_subtraction:
            pulse = self.curr_leak_amp(I, dur1, dur2, delay, dt)
            current = current - pulse

        if chord_conductance:
            current = current/(self.vclamp.amp2 - self.E_rev)
        return current

    @staticmethod
    def normalize_to_one(current, normalization="to_one"):
        """
        normalization: "to_one" or "save_sign"
        """
        values = np.array(list(current.values()))
        factor = max(abs(values))
        new_current = {}
        if normalization == "save_sign":
            for key in current.keys():
                new_current[key] = current[key]/factor
        else:
            for key in current.keys():
                new_current[key] = abs(current[key])/factor
        return new_current

    @staticmethod
    def get_max_of_dict(current, ion_name, chord_conductance):
        new_current = {}
        for key in current.keys():
            #inward currents are negative
            if ion_name in ["Ca", "ca", "Ba", "ba", "na"]:
                if not chord_conductance:
                    new_current[key] = current[key].min()
                else:
                    new_current[key] = current[key].max()
            elif ion_name == "k":
                new_current[key] = current[key].max()
            else:
                new_current[key] = abs(current[key]).max()
        return new_current


class ModelPatchNernst(ModelPatchWithChannel):    
    def __init__(self, path_to_mods, channel_name, ion_name,
                 external_conc=None, gbar_name="gbar", temp=22, recompile=True,
                 liquid_junction_pot=0, cvode=True, R_m=20000, v_rest=-65,
                 E_rev=None, gbar_value=None):
        super(ModelPatchNernst, self).__init__(path_to_mods, channel_name,
                                               ion_name,
                                               external_conc=external_conc,
                                               gbar_name=gbar_name, temp=temp,
                                               recompile=recompile,
                                               liquid_junction_pot=liquid_junction_pot,
                                               cvode=cvode,
                                               R_m=20000, v_rest=v_rest,
                                               gbar_value=gbar_value)
        if ion_name.lower() == "nonspecific":
            if isinstance(E_rev, int) or isinstance(E_rev, float):
                self.E_rev = E_rev
            else:
                raise SystemExit('Unable to proceed, if E_rev is unknown.')
        else:
            E_rev_name = self._find_E_rev_name()
            self.E_rev = self.calc_E_rev(E_rev)
            setattr(self.patch,  E_rev_name, self.E_rev)

    
    @property
    def cai(self):
        return self._cai
   
    @cai.setter
    def cai(self, value):
        self._cai = value
        if self.ion_name.lower() == "ca" and self._external_conc is not None:
            self.E_rev = self.calc_E_rev()

    @property
    def Cai(self):
        return self._cai

    @cai.setter
    def Cai(self, value):
        self._cai = value
        if self.ion_name.lower() == "ca" and self._external_conc is not None:
            self.E_rev = self.calc_E_rev()

    @property
    def external_conc(self):
        return self._external_conc

    @external_conc.setter
    def external_conc(self, value):
        self._external_conc = value
        self.E_rev = self.calc_E_rev()
    @property
    def ki(self):
        return self._ki

    @ki.setter
    def ki(self, value):
        self._ki = value
        if self.ion_name == "k" and self._external_conc is not None:
            self.E_rev = self.calc_E_rev()

    @property
    def nai(self):
        return self._nai

    @nai.setter
    def nai(self, value):
        self._nai = value
        if self.ion_name == "na" and self._external_conc is not None:
            self.E_rev = self.calc_E_rev()


    
class ModelPatchCaConcentration(ModelPatchWithChannel):
    def __init__(self, path_to_mods, channel_name, ion_name,
                 Ca_ext, gbar_name="gbar", temp=22, recompile=True,
                 liquid_junction_pot=0, cvode=True, R_m=20000, v_rest=-65,
                 gbar_value=None, t_decay=20, buffer_capacity=18,
                 membrane_shell_width=membrane_shell_width):
        """
        Model class for testing calcium channels with 
        """
        super(ModelPatchCaConcentration, self).__init__(path_to_mods,
                                                      channel_name,
                                                      ion_name,
                                                      external_conc=None,
                                                      gbar_name=gbar_name,
                                                      temp=temp,
                                                      recompile=recompile,
                                                      liquid_junction_pot=liquid_junction_pot,
                                                      cvode=cvode,
                                                      R_m=20000, v_rest=v_rest,
                                                      gbar_value=gbar_value)
        self.E_rev = None
        self.t_decay = t_decay
        self.Kb = buffer_capacity
        self._Ca_ext = Ca_ext
        self.memb_shell_width = membrane_shell_width
        self.geom = rxd.Shell(1- self.memb_shell_width, 1)
        self.memb_shell = rxd.Region(self.patch,
                                     geometry=self.geom,
                                     nrn_region="i",
                                     name="membrane_shell")
        old_saf = self.geom.surface_areas1d
        self.geom.surface_areas1d = lambda sec: old_saf(sec)/self.Kb
        if self.ion_name == "ca":
            self.ca = rxd.Species(self.memb_shell, d=0.2,
                                  name='ca', charge=2,
                                  initial=self._cai,
                                  atolscale=1e-9)
            h.cao0_ca_ion = self._Ca_ext
            self.E_Ca_rev = self.calc_E_rev(external=self.Ca_ext)
        elif self.ion_name == "Ca":
            self.ca = rxd.Species(self.memb_shell, d=0.2,
                                  name='Ca', charge=2,
                                  initial=self._cai,
                                  atolscale=1e-9)
            h.Cao0_Ca_ion = self._Ca_ext
            self.E_Ca_rev = self.calc_E_rev(external=self._Ca_ext)
        elif self.ion_name == "Ba":
            self.ca = rxd.Species(self.memb_shell, d=0.2,
                                  name='Ca', charge=2,
                                  initial=0,
                                  atolscale=1e-9)
            self._cai = 0
            self.E_Ca_rev = None
            h.Cao0_Ca_ion = self._Ca_ext
            chan = self.patch.psection()["density_mechs"][channel_name]
            for i, seg in enumerate(self.patch):
                from_mech = getattr(seg, channel_name)
                gbar_val = 2*chan[gbar_name][i]
                setattr(from_mech, gbar_name, gbar_val)
        elif self.ion_name == "ba":
            self.ca = rxd.Species(self.memb_shell, d=0.2,
                                  name='ca', charge=2,
                                  initial=0,
                                  atolscale=1e-9)
            self._cai = 0
            self.E_Ca_rev = None
            h.cao0_ca_ion = self._Ca_ext
            chan = self.patch.psection()["density_mechs"][channel_name]
            for i, seg in enumerate(self.patch):
                from_mech = getattr(seg, channel_name)
                gbar_val = 2*chan[gbar_name][i]
                setattr(from_mech, gbar_name, gbar_val)
        else:
            raise SystemExit("Unknown ion %s. I only know Ca, ca, Ba and ba"
                             % self.ion_name)
        self.decay_eq = (self._cai - self.ca)/self.t_decay
        self.ca_decay = rxd.Rate(self.ca, self.decay_eq)

    @property
    def cai(self):
        if self.ion_name == "ca":
            return self.patch.cainf_cad
        elif self.ion_name == "Ca":
            return self.patch.cainf_Cad
        elif self.ion_name == "ba":
            return 0
        elif self.ion_name == "Ba":
            return 0

    @cai.setter
    def cai(self, value):
        self._cai = value
        if self.ion_name.lower() == "ca":
            self.patch.cainf_cad = self._cai
            self.E_Ca_rev = self.calc_E_rev(external=self._Ca_ext)
        elif self.ion_name == "Ca":
            self.patch.cainf_Cad = self._cai
            self.E_Ca_rev = self.calc_E_rev(external=self._Ca_ext)
                             
    @property
    def Cai(self):
        if self.ion_name == "ca":
            return self.patch.cainf_cad
        elif self.ion_name == "Ca":
            return self.patch.cainf_cad
        else:
            return 0

    @cai.setter
    def Cai(self, value):
        self._cai = value
        if self.ion_name.lower() == "ca":
            self.patch.cainf_cad = self._cai
            self.E_Ca_rev = self.calc_E_rev(external=self._Ca_ext)
        elif self.ion_name == "Ca":
            self.patch.cainf_Cad = self._cai
            self.E_Ca_rev = self.calc_E_rev(external=self._Ca_ext)

    @property
    def Ca_ext(self):
        return  self._Ca_ext

    @Ca_ext.setter
    def Ca_ext(self, value):
        self._Ca_ext = value
        if self.ca_ion == "ca":
            h.cao0_ca_ion = self._Ca_ext
            self.E_Ca_rev = self.calc_E_rev(external=self._Ca_ext)
        elif self.ca_ion == "Ca":
            h.Cao0_Ca_ion = self._Ca_ext
            self.E_Ca_rev = self.calc_E_rev(external=self._Ca_ext)
        elif self.ca_ion == "Ba":
            h.Cao0_Ca_ion = self._Ca_ext
        elif self.ca_ion == "ba":
            h.cao0_ca_ion = self._Ca_ext
        else:
            raise SystemExit("Unknown ion %s. I only know Ca and ca"
                             % self.ion_name )

class WholeCellAttributes:
    
    def _set_g_pas(self, R_in, sec_list):
        self._R_in = R_in
        area = 0
        for sec in sec_list:
            for seg in sec:
                area += seg.area()*1e-4 #  in cm2
        for sec in sec_list:
            sec.g_pas = 1/(self._R_in*area)  # in mho/cm2
            
    @property
    def R_in(self):
        return self._R_in

    @R_in.setter
    def R_in(self, value):
        self._set_g_pas(value)


class ModelWholeCellPatchNernst(ModelPatchNernst, WholeCellAttributes):
    """
    R_in -- in ohms
    """
    def __init__(self, path_to_mods, channel_name, ion_name,
                 external_conc=None, gbar_name="gbar",
                 temp=22, recompile=True,
                 liquid_junction_pot=0, cvode=True,  R_in=200e6,
                 v_rest=-65, E_rev=None, gbar_value=None):

        super(ModelWholeCellPatchNernst, self).__init__(path_to_mods,
                                                        channel_name,
                                                        ion_name,
                                                        external_conc,
                                                        gbar_name, temp,
                                                        recompile,
                                                        liquid_junction_pot,
                                                        cvode, 20000, v_rest,
                                                        E_rev, gbar_value)
        self._L = 10
        self._diam = 10
        self.patch.L = self._L
        self.patch.Ra = 100
        self.patch.diam = self._diam
        self._set_g_pas(R_in, [self.patch])

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


class ModelWholeCellPatchConcentration(ModelPatchCaConcentration,
                                       WholeCellAttributes):
    """
    R_in -- in ohms
    """
    def __init__(self, path_to_mods, channel_name, ion_name,
                 external_conc=None, gbar_name="gbar",
                 temp=22, recompile=True,
                 liquid_junction_pot=0, cvode=True,  R_in=200e6,
                 v_rest=-65, gbar_value=None):

        super(ModelWholeCellPatchConcentration, self).__init__(path_to_mods,
                                                               channel_name,
                                                               ion_name,
                                                               external_conc,
                                                               gbar_name, temp,
                                                               recompile,
                                                               liquid_junction_pot,
                                                               cvode, 20000,
                                                               v_rest,
                                                               gbar_value)
        self._L = 10
        self._diam = 10
        self.patch.L = self._L
        self.patch.Ra = 100
        self.patch.diam = self._diam
        self._set_g_pas(R_in, [self.patch])

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


