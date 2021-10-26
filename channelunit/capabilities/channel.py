import sciunit


class NModlChannel(sciunit.Capability):
    """This class provides capabilities for a ion channel implemented
    using nmodl. For now it is calculation of the activation and
    inactivation steady state curves.

    """


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
