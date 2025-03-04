import os
from subprocess import run

import numpy as np
from scipy.signal import bessel, lfilter

import sciunit
from neuron import h
from neuron import rxd
import neuron
from channelunit.capabilities import NModlChannel


loc = os.path.dirname(os.path.abspath(__file__))
mechanisms_path = os.path.join(loc, 'mechanisms')



N = 4
DT = 5e-5  #  Magee digitizes at 20kHz


F = 96485.33212  # C mol^-1
R = 8.314462618  # J mol^-1 K^-1

class MembranePatch(sciunit.Model):
    def __init__(self, temp, Rm, cm,
                 v_rest, ljp, cvode, sim_dt):
        h.load_file("stdrun.hoc")
        self.dt = DT
        self.compile_and_add(mechanisms_path, True)
        self.junction = ljp
        self.patch = h.Section(name="patch")
        self.patch.L = 1
        self.patch.Ra = 100
        self.patch.diam = 1
        self.patch.insert("pas")
        self.patch.e_pas = v_rest
        self.patch.g_pas = 1/Rm
        self.patch.cm = cm # uF/cm2
        self.temperature = temp
        self.v_low = 50
        self.vclamp = h.SEClampOLS(self.patch(0.5))
        if cvode is True:
            self.cvode = True
            h.CVode()
            h.CVode().atol(1e-7)
        else:
            self.sim_dt = sim_dt
            h.dt = self.sim_dt
            self.cvode = False
        self.f_b, self.f_a = bessel(8, 100, btype='low',
                                    analog=False, norm='mag', fs=1/DT)
        
    def compile_and_add(self, path, recompile):
        working_dir = os.getcwd()
        os.chdir(path)
        if recompile:
            p = run('nrnivmodl')
        neuron.load_mechanisms(path)
        os.chdir(working_dir)

    def set_vclamp(self, dur1, v1, dur2, v2, leak_subtraction):
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
            return dur1+dur2

        pulse_amp = self.vclamp.amp2 - self.vclamp.amp1
        v_sub = self.vclamp.amp1 - pulse_amp/N - self.v_low 
        v_pulse = v_sub + pulse_amp/N 

        self.vclamp.dur3 = dur2
        self.vclamp.amp3 = self.vclamp.amp1
        self.vclamp.amp4 = v_sub
        self.vclamp.dur4 = dur2
        # 1st pulse
        self.vclamp.amp5 = v_pulse
        self.vclamp.dur5 = dur2
        self.vclamp.amp6 = v_sub
        self.vclamp.dur6 = dur2
        #2nd pulse
        self.vclamp.amp7 = v_pulse
        self.vclamp.dur7 = dur2
        self.vclamp.amp8 = v_sub
        self.vclamp.dur8 = dur2
        #3rd pulse
        self.vclamp.amp9 = v_pulse
        self.vclamp.dur9 = dur2
        self.vclamp.amp10 = v_sub
        self.vclamp.dur10 = dur2
        #4th pulse
        self.vclamp.amp11 = v_pulse
        self.vclamp.dur11 = dur2
        self.vclamp.amp12 = v_sub
        self.vclamp.dur12 = dur2
        return dur1 + 11*dur2

    def run(self, t_stop, dt=DT):
        current = h.Vector()
        current.record(self.vclamp._ref_i, dt)
        h.finitialize(self.patch.e_pas)
        if self.cvode:
            h.CVode().re_init()
        else:
            h.fcurrent()
        h.frecord_init()
        h.continuerun(t_stop)
        return current.as_numpy()

    @classmethod
    def curr_stim_response(self, I, dur1, dur2, dt,
                           leak_subtraction):
        beg_stim = int(np.round(dur1/dt))
        end_stim = int(np.round((dur1+dur2)/dt))
        current = I.copy()
        if leak_subtraction:
            beg = int(np.round(abs(dur1 - dur2)/dt))
            current[beg_stim: end_stim] -= I[beg: beg_stim]
        return current

    @classmethod
    def curr_leak_amp(self, I, dur1, dur2, dt):
        t_start =  dur1 + 2*dur2 
        length = int(dur2/dt)
        pulse = np.zeros((length,))
        for i in range(N):
            bas_start = int(np.round(t_start/dt))
            pulse_start = int(np.round((t_start+dur2)/dt))
            basal = I[bas_start:pulse_start].copy()
            pulse_end = int(np.round((t_start+2*dur2)/dt))
            pulse += (I[pulse_start:pulse_end].copy() - basal)
            t_start += 2*dur2
        return pulse.copy()

    @property
    def cm(self):
        return self.patch.cm

    @cm.setter
    def cm(self, value):
        self.patch.cm = value

    @property
    def Rm(self):
        return 1/self.patch.g_pas

    @Rm.setter
    def Rm(self, value):
        assert value > 0
        self.patch.g_pas = 1/value


class ModelPatch(MembranePatch, NModlChannel):
    _nai = 10  # mM
    _ki = 140  # mM
    _cai = 50e-6  # mM
    def __init__(self, path_to_mods: str, channel_names: list, ion_names: list,
                 external_conc: dict, internal_conc: dict, E_rev: dict,
                 gbar_names: dict, gbar_values: dict,
                 temp, recompile: bool, ljp, cvode: bool,
                 Rm, cm, v_rest, directory,
                 sim_dt):
        """
        ion_name: list
            most common ions are: na (sodium), k, ca. 
            if you specify nonspecific, you need to provide a value of E_rev
            It is important to pay attention whether the ion variable is 
            specified with a lowercase or an uppercase letter, because 
            the name of the reversal potential variable is constructed
            based on the ion name (na -> ena, k -> ek, ca -> eca).
        Rm: float 
            mebrane resistivity (in ohm*cm^2)
        """
        if not directory:
            directory = "validation_results"
        self.channel_names = []
        self.mod_path = path_to_mods
        self.compile_and_add(self.mod_path, recompile)
        super(ModelPatch, self).__init__(temp=temp, 
                                         Rm=Rm, cm=cm,
                                         v_rest=v_rest,
                                         ljp=ljp,
                                         cvode=cvode,
                                         sim_dt=sim_dt,)

        self.channels = []
        self.gbar_names = {}
        self.external_conc = {"Ca": 0}
        self.E_rev = {}
        self.base_directory = directory
        self.ion_names = ion_names
        #set up channel conductance/permeability in case it is 0
        for channel in channel_names:
            if channel not in gbar_names:
                gbar_name = "gbar"
            else:
                gbar_name = gbar_names[channel]
            if channel not in gbar_values:
                gbar_value = 0.001
            else:
                gbar_value = gbar_values[channel]

            self.add_channel(channel, gbar_name, gbar_value)

        for ion in ion_names:
            if ion in external_conc:
                e_conc = external_conc[ion]
                if ion == "ca" or ion == "ba":
                    self.external_conc["Ca"] = e_conc
                else:
                    self.external_conc[ion] = e_conc
            else:
                e_conc = None
            if ion in internal_conc:
                if ion.lower() == "ca" or ion.lower() == "ba":
                    self._cai = internal_conc[ion]
                elif ion == "k":
                    self._ki = internal_conc[ion]
                elif ion == "na":
                    self._nai = internal_conc[ion]

            if ion in E_rev:
                e_rev = E_rev[ion]
            else:
                e_rev = None
                
            self.E_rev[ion] = self.calc_E_rev(ion, e_rev, e_conc)
            if ion == "k":
                self.patch.ek = self.E_rev[ion]
            elif ion == "na":
                self.patch.ena = self.E_rev[ion]
            elif ion == "ca":
                self.patch.eca = self.E_rev[ion]
        self.ca = None
        
    def add_channel(self, channel_name, gbar_name, gbar_value):
        self.channel_names.append(channel_name)
        self.patch.insert(channel_name)

        chan = self.patch.psection()["density_mechs"][channel_name]
        self.channels.append(chan)
        if gbar_name not in chan.keys():
            raise SystemExit('Unable to proceed, unknown %s conductance (gbar)'
                             % channel_name)
        self.gbar_names[channel_name] = gbar_name
        if chan[gbar_name][0] == 0:
            for seg in self.patch:
                from_mech = getattr(seg, channel_name)
                setattr(from_mech, gbar_name, gbar_value)


    def get_gbar(self, channel_name):
        mech = self.patch.psection()["density_mechs"][channel_name]
        values = mech[self.gbar_names[channel_name]]
        if len(values) == 1:
            return values[0]
        elif len(set(values)) == 1:
            return values[0]
        return values


    def set_gbar(self, channel_name, value):
        if not isinstance(value, list):
            value = [value]
        if len(value) == self.patch.nseg:
            for i, seg in enumerate(self.patch):
                from_mech = getattr(seg, channel_name)
                setattr(from_mech, self.gbar_names[channel_name], value[i])
        elif len(value) == 1:
            for i, seg in enumerate(self.patch):
                from_mech = getattr(seg, channel_name)
                setattr(from_mech, self.gbar_names[channel_name], value[0])
       
    def calc_E_rev(self, ion_name, E_rev=None, external=None):
        if ion_name.lower() == "k":
            internal = self._ki
            valence = 1
        elif ion_name.lower() == "na":
            internal = self._nai
            valence = 1
        elif ion_name.lower() == "ca":
            internal = self._cai
            valence = 2
        elif ion_name.lower() == "ba":
            return
        elif E_rev is None:
            raise SystemExit("Unknown ion type %s. Only now na, k, ca"
                              % ion_name)
        elif external is not None:
            raise SystemExit("Unknown ion type %s. Only now na, k, ca"
                             % ion_name)

        if external is None:
            if E_rev is None:
                name = "e%s" % ion_name
                return self.patch.psection()["ions"][ion_name][name][0]
            else:
                return E_rev
        elif not isinstance(external, int) and not isinstance(external, float):
            if E_rev is None:
                name = "e%s" % ion_name
                return self.patch.psection()["ions"][ion_name][name][0]
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
        for channel_name in self.channel_names:
            fname = "%s_%s" % (fname, channel_name)
        for ion_name in self.external_conc.keys():
            if self.external_conc[ion_name]:
                fname = "%s_%s_ext_%f_mM" % (fname, ion_name,
                                                self.external_conc[ion_name])
        if "ca" in self.ion_names:
             fname = "%s_int_Ca_%f_mM" % (fname, self._cai)
        if chord_conductance:
            for ion_name in self.E_rev:
                if self.E_rev[ion_name] is not None:
                    fname = "%s_E_%s_rev_%f_mV" % (fname, ion_name,
                                                      self.E_rev[ion_name])
      
        fname = "%s_from_%4.2f_mV_to_%f_mV" % (fname, stim_beg, stim_end)
        
        if self.cvode:
            fname = "%s_cvode" % fname
        if ca_conc:
            fname = "%s_Ca_conc" % fname
        return fname

    def _record_current(self, electrode_current, chord_conductance):
        current = h.Vector()
        if electrode_current:
            current.record(self.vclamp._ref_i, self.dt)
        else:
            if len(self.ion_names) > 1:
                current.record(self.vclamp._ref_i, self.dt)
                chord_conductance = False
            else:
                if self.ion_names[0] == "na":
                    current.record(self.patch(0.5)._ref_ina, self.dt)
                elif self.ion_names[0] == "k":
                    current.record(self.patch(0.5)._ref_ik, self.dt)
                elif self.ion_names[0] == "ca":
                    current.record(self.patch(0.5)._ref_ica, self.dt)
                elif self.ion_names[0] == "ba":
                    current.record(self.patch(0.5)._ref_ica, self.dt)
                else:
                    current.record(self.vclamp._ref_i, self.dt)
        return current, chord_conductance

    def get_activation_traces(self, stimulation_levels: list,
                              v_hold: float, t_stop:float,
                              chord_conductance=False,
                              electrode_current=True,
                              save_traces=True, save_ca=False):
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
        out_I = []
        delay = t_stop
        shift = 0 #(self.patch.cm/self.patch.g_pas)*1e-3
        if not electrode_current:
            leak_subtraction = False
        else:
            leak_subtraction = True
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
            cal_ref = self.ca[self.memb_shell].nodes[0]._ref_concentration
            calcium.record(cal_ref,
                           self.dt)
        else:
            ca_fname = ""
        h.celsius = self.temperature
        current_vals = {}
        current,\
            chord_conductance = self._record_current(electrode_current,
                                                     chord_conductance)
        time = h.Vector()
        time.record(h._ref_t, self.dt)
        voltage = h.Vector()
        voltage.record(self.patch(0.5)._ref_v, self.dt)
        stim_start = int(delay/self.dt)
        out_v = []
        beg = int(np.round((shift+delay)/self.dt))
        end = int(np.round((delay+t_stop)/self.dt))
        for i, level in enumerate(stimulation_levels):
            stim_stop = self.set_vclamp(delay, v_hold, t_stop, level,
                                        leak_subtraction)
            h.finitialize(v_hold)
            if self.cvode:
                h.CVode().re_init()
            else:
                h.fcurrent()
            h.frecord_init()

            h.continuerun(stim_stop)

            I = current.as_numpy()
            if not i:
                out_I.append(time.as_numpy())
            out_I.append(current.as_numpy().copy())
            out = self.extract_current(I, chord_conductance,
                                       leak_subtraction,
                                       delay+shift,t_stop, self.dt,
                                       level)

            if save_ca:
                if not i:
                    save_time = time.as_numpy().copy()
                    calcium_vals.append(save_time[beg: end])
  
                calcium_vals.append(calcium.as_numpy().copy()[beg: end])
            if save_traces:
                if not i:
                    save_time = time.as_numpy().copy()
                    output.append(save_time[beg: end].copy())
                output.append(out[beg: end].copy())
                
                header += ";%4.2f" % level
            current_vals[level] = out[beg: end].copy()

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

    def get_activation_SS(self, stimulation_levels: list,
                          v_hold: float, t_stop:float,
                          power: int, t_mes, chord_conductance,
                          electrode_current,
                          normalization="to_one",
                          save_traces=True, save_ca=False):
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
        if len(self.channel_names) > 1:
            electrode_current = True
            chord_conductance = False
        currents = self.get_activation_traces(stimulation_levels,
                                              v_hold, t_stop,
                                              chord_conductance,
                                              electrode_current,
                                              save_traces=save_traces,
                                              save_ca=save_ca)
        
        max_current = self.get_max_of_dict(currents, t_mes, self.dt)
        result = self.normalize_to_one(max_current, normalization)
        if power != 1:
            for key in result.keys():
                result[key] = result[key]**(1/power)
        return result

    def get_inactivation_traces(self, stimulation_levels: list,
                                v_test: float, t_test:float,
                                chord_conductance,
                                electrode_current,
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
        delay = t_test
        shift = 0 # (self.patch.cm/self.patch.g_pas)*1e-3

        if not electrode_current:
            leak_subtraction = False
        else:
            leak_subtraction = True
        
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
        stim_start = int(delay/self.dt)
        current_values = {}
        current,\
            chord_conductance = self._record_current(electrode_current,
                                                     chord_conductance)
        time = h.Vector()
        time.record(h._ref_t, self.dt)
        voltage = h.Vector()
        voltage.record(self.patch(0.5)._ref_v, self.dt)

        for i, v_hold in enumerate(stimulation_levels):
            t_stop = self.set_vclamp(delay, v_hold, t_test, v_test,
                                     leak_subtraction)
            h.finitialize(v_hold)
            if self.cvode:
                h.CVode().re_init()
            else:
                h.fcurrent()
            h.frecord_init()
            h.continuerun(t_stop)
            I = current.as_numpy()
            out = self.extract_current(I,
                                       chord_conductance,
                                       leak_subtraction, delay,
                                       t_test, self.dt,
                                       v_hold)

            beg = int(np.round(delay/self.dt))
            end = int(np.round((delay+t_test)/self.dt))
            if save_ca:
                if not i:
                    save_time = time.as_numpy().copy()
                    calcium_vals.append(save_time[beg: end])
                calcium_vals.append(calcium.as_numpy().copy()[beg: end])
            if save_traces:
                if not i:
                    save_time = time.as_numpy().copy()
                    output.append(save_time[beg: end])
                output.append(out[beg: end].copy())
                header += ";%4.2f" % v_hold
            current_values[v_hold] = out.copy()
        if save_traces:
            path = os.path.join(self.base_directory, "data")
            if not os.path.exists(path):
                makedirs(path)
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
        
    def get_inactivation_SS(self, stimulation_levels: list,
                            v_test: float, t_test:float,
                            power: int, t_mes,
                            chord_conductance=False,
                            electrode_current=True,
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
        if len(self.channel_names) > 1:
            electrode_current = True
            chord_conductance = False
        currents = self.get_inactivation_traces(stimulation_levels,
                                                v_test, t_test,
                                                chord_conductance,
                                                electrode_current,
                                                save_traces=save_traces,
                                                save_ca=save_ca)
        max_current = self.get_max_of_dict(currents, t_mes, self.dt)
        result = self.normalize_to_one(max_current, normalization)
        if power != 1:
            for key in result.keys():
                result[key] = result[key]**(1/power)
        return result
                
    def extract_current(self, I, chord_conductance, leak_subtraction, dur1,
                        dur2, dt, v=None, filtering=True):
        if v is None:
            v= self.patch.e_pas
        if filtering:
            filter_current = lfilter(self.f_b, self.f_a, I)
        else:
            filter_current = I.copy()
        current = self.curr_stim_response(filter_current, dur1, dur2, dt,
                                          leak_subtraction)
        #either step injection or the short pulse
        if leak_subtraction:
            beg = int(dur1/dt)
            end = int((dur1+dur2)/dt)
            pulse = self.curr_leak_amp(filter_current, dur1, dur2, dt)
            current[beg:end] = current[beg:end] - pulse
        if chord_conductance:
            current = current/(v - self.E_rev[self.ion_names[0]])
        return current

    @staticmethod
    def normalize_to_one(current, normalization="save_sign"):
        """
        normalization: "to_one" or "save_sign"
        """
        values = np.array(list(current.values()))
        factor = max(abs(values))
        new_current = {}
        if factor < 1e-2:
            for key in current.keys():
                new_current[key] = 0
        if normalization == "save_sign":
            for key in current.keys():
                new_current[key] = current[key]/factor
        else:
            for key in current.keys():
                new_current[key] = abs(current[key])/factor
        return new_current

    @staticmethod
    def get_max_of_dict(current, t_mes, dt):
        new_current = {}
        for key in current.keys():
            if t_mes is None:
                #inward currents are negative
                if max(current[key]) < max(abs(current[key])):
                    new_current[key] = min(current[key])                
                else:
                    new_current[key] = max(current[key])
            else:
                 new_current[key] = current[key][int(t_mes/dt)]
        return new_current

    
    @property
    def cai(self):
        return self._cai
   
    @cai.setter
    def cai(self, value):
        self._cai = value
        if self.external_conc["ca"] is not None:
            self.E_rev["ca"] = self.calc_E_rev("ca", None,
                                               self.external_conc["Ca"])

    @property
    def Cai(self):
        return self._cai

    @cai.setter
    def Cai(self, value):
        self._cai = value
       
        if self.external_conc["ca"] is not None:
            self.E_rev["ca"] = self.calc_E_rev("ca", None,
                                               self.external_conc["Ca"])
            self.patch.eca = self.E_rev["ca"]

    def get_external_conc(self, ion):
        return self.external_conc[ion]

    def set_external_conc(self, ion, value):
        self.external_conc[ion] = value
        self.E_rev[ion] = self.calc_E_rev(ion, None,
                                          self.external_conc[ion])
        if ion == "k":
            self.patch.ek = self.E_rev[ion]
            
        elif ion == "na":
            self.patch.ena = self.E_rev[ion]
        elif ion == "ca":
            self.patch.eca = self.E_rev[ion]

    @property
    def ki(self):
        return self._ki

    @ki.setter
    def ki(self, value):
        self._ki = value
        if self.external_conc["k"] is not None:
            self.E_rev["k"] = self.calc_E_rev("k", None,
                                              self.external_conc["k"])
            self.patch.ek = self.E_rev["k"]

    @property
    def nai(self):
        return self._nai

    @nai.setter
    def nai(self, value):
        self._nai = value
        if self.external_conc["na"] is not None:
            self.E_rev["na"] = self.calc_E_rev("na", None,
                                               self.external_conc["na"])
            self.patch.ena = self.E_rev["na"]


class ModelPatchCa(ModelPatch):
    def __init__(self, path_to_mods: str, channel_names: list, ion_names: list,
                 external_conc: dict, internal_conc: dict, E_rev: dict,
                 gbar_names: dict, gbar_values: dict,
                 temp, recompile: bool, ljp, cvode: bool,
                 Rm, cm, v_rest, directory,
                 sim_dt, t_decay, L, diam, Ra,
                 buffer_capacity,
                 membrane_shell_width):

        if not directory:
            directory = "validation_results"
        super(ModelPatchCa, self).__init__(path_to_mods, channel_names,
                                           ion_names,
                                           external_conc=external_conc,
                                           internal_conc=internal_conc,
                                           E_rev=E_rev,
                                           gbar_names=gbar_names,
                                           gbar_values=gbar_values,
                                           temp=temp, recompile=recompile,
                                           ljp=ljp, cvode=cvode, Rm=Rm,
                                           cm=cm, v_rest=v_rest,
                                           directory=directory, sim_dt=sim_dt)

        self.t_decay = t_decay
        self.Kb = buffer_capacity
        self.memb_shell_width = membrane_shell_width
        self.geom = rxd.Shell(1- self.memb_shell_width, 1)
        self.memb_shell = rxd.Region(self.patch,
                                     geometry=self.geom,
                                     nrn_region="i",
                                     name="membrane_shell")
        old_saf = self.geom.surface_areas1d
        self.geom.surface_areas1d = lambda sec: old_saf(sec)/self.Kb
        if "ca" in self.ion_names:
            self.ca = rxd.Species(self.memb_shell, d=0.2,
                                  name='ca', charge=2,
                                  initial=self._cai,
                                  atolscale=1e-9)
            h.cao0_ca_ion = self.external_conc["Ca"]

        elif "ba" in self.ion_names:
            self.ca = rxd.Species(self.memb_shell, d=0.2,
                                  name='ca', charge=2,
                                  initial=0,
                                  atolscale=1e-9)
            if "ba" not in internal_conc:
                self._cai = 0
            self.E_rev["ca"] = None
            h.cao0_ca_ion = self.external_conc["Ca"]
            for channel_name in self.channel_names:
                chan = self.patch.psection()["density_mechs"][channel_name]
                for i, seg in enumerate(self.patch):
                    from_mech = getattr(seg, channel_name)
                    gbar_val = 2*chan[self.gbar_names[channel_name]][i]
                    setattr(from_mech, self.gbar_names[channel_name], gbar_val)
        
        self.decay_eq = (self._cai - self.ca)/self.t_decay
        self.ca_decay = rxd.Rate(self.ca, self.decay_eq)

    @property
    def cai(self):
        if "ca" in self.ion_names:
            return self.patch._cai
        elif "ba" in self.ion_names:
            return 0

    @cai.setter
    def cai(self, value):
        self._cai = value
        if value > 0:
            if "ca" in self.ion_names:
                self.patch.eca = self.calc_E_rev("ca",
                                                 external=self.external_conc["Ca"])
            
    @property
    def Ca_ext(self):
        return  self.external_conc["Ca"]

    @Ca_ext.setter
    def Ca_ext(self, value):
        self._Ca_ext = value
        if self.ca_ion == "ca":
            h.cao0_ca_ion = self._Ca_ext
            self.patch.eca = self.calc_E_rev("ca",
                                            external=self.external_conc["Ca"])
        elif self.ca_ion == "ba":
            h.cao0_ca_ion = self.external_conc["Ca"]
        else:
            raise SystemExit("Unknown ion %s. I only know ca"
                             % self.ion_name )
        
class WholeCellAttributes:

    def area(self, sec_list):
        area = 0
        for sec in sec_list:
            for seg in sec:
                area += seg.area()
        return area

    def _set_g_pas(self, Rin, sec_list):
        self._Rin = Rin
        area = self.area(sec_list)*1e-8 # in cm2
        for sec in sec_list:
            sec.g_pas = 1/(self._Rin*area)  # in mho/cm2

    def set_cap(self, cap, sec_list):
        self._cap = cap
        area = self.area(sec_list)*1e-8 #in cm2        
        for sec in sec_list:
            sec.cm = self.cap/area*1e5
        # in uF/cm2, starting from nF, 1E-9F/1E-8cm2=1E-1
        # F/cm2=1e5 1E-6F/cm2=10 uF/cm2
        
    @property
    def Rin(self):
        return self._Rin

    @Rin.setter
    def Rin(self, value):
        self._set_g_pas(value, [self.patch])

    @property
    def cap(self):
        return self._cap
    
    @cap.setter
    def cap(self, value):
        self.set_cap(value, [self.patch])

    @property
    def L(self):
        return self.patch.L

    @L.setter
    def L(self, value):
        self._L = value
        self.patch.L = self._L
        self._set_g_pas(self.Rin, [self.patch])

    @property
    def diam(self):
        return self.patch.diam

    @diam.setter
    def diam(self, value):
        self._diam = value
        self.patch.diam = self._diam
        self._set_g_pas(self.Rin,  [self.patch])


