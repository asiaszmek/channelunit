import sciunit


class NModlChannel(sciunit.Capability):
    """
    Th
    """

    def get_E_rev(self, soma):
        ions = soma.psection()["ions"]
        if not len(ions):
            return None
        if len(ions) > 1:
            return None
        return "e%s" % ions[0]

    def get_activation_steady_state(self, soma_obj, clamp_obj,
                                    stimulation_levels: list,
                                    v_init: float, t_stop:float,
                                    chord_conductance=False,
                                    E_rev=None,
                                    duration=400):
        """
        channel should be a neuron density mechanism, clamp a SEClamp neuron 
        object
        stimulation dict {amplitude: stim_duration}
        """
        dt = 0.01
        max_current = {}
        current = h.Vector().record(soma_obj(0.5)._ref_i_membrane, dt)
        time = h.Vector().record(h._ref_t, dt)
        vm = h.Vector().record(soma_obj(0.5)._ref_v, dt)
        delay = 200
        if chord_conductance and E_rev is None:
            E_rev = self.get_E_rev(soma_obj)
            if E_rev is None:
                raise SystemExit('Unable to calculate chord conductance, if E_rev is unknown.')

        for level in stimulation_levels:
            clampobj.dur1 = delay
            clampobj.amp1 = v_init
            clampobj.dur2 = duration
            clampobj.amp2 = level
            h.dt = 0.001
            h.tstop = t_stop
            h.run(t_stop)
            if chord_conductance:
                out = current.as_numpy()/(vm - E_rev)
            else:
                out = abs(current.as_numpy())
            max_current[level] = max(out)
        return max_current

    def get_inactivation_steady_state(self, soma_obj, clamp_obj,
                                      stimulation_levels: list,
                                      v_test: float, t_test:float,
                                      chord_conductance=False):
        dt = 0.01
        delay = 200
        t_stop = delay +  t_test
        max_current = {}
        current = h.Vector().record(soma_obj(0.5)._ref_i_membrane, dt)
        time = h.Vector().record(h._ref_t, dt)
        vm = h.Vector().record(soma_obj(0.5)._ref_v, dt)
        if chord_conductance and E_rev is None:
            E_rev = self.get_E_rev(soma_obj)
            if E_rev is None:
                raise SystemExit('Unable to calculate chord conductance, if E_rev is unknown.')
        for level in stimulation_levels:
            clampobj.dur1 = delay
            clampobj.amp1 = level
            clampobj.dur2 = t_test
            clampobj.amp2 = v_test
            h.dt = 0.001
            h.tstop = t_stop
            h.run(t_stop)
            if chord_conductance:
                out = current.as_numpy()/(vm - E_rev)
            else:
                out = abs(current.as_numpy())
            max_current[level] = max(out)
        return max_current
