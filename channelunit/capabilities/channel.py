import sciunit

class NModlChannel(sciunit.Capability):
    """
    Th
    """
    def get_activation_steady_state(self, soma_obj, clamp_obj,
                                    stimulation_levels: list,
                                    v_init: float, t_stop:float,
                                    chord_conductance=False, E_rev=0,
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
        for level in stimulation_levels:
            clampobj.dur1 = delay
            clampobj.amp1 = v_init
            clampobj.dur2 = duration
            clampobj.amp2 = level
            h.dt = 0.001
            h.tstop = t_stop
            h.run(t_stop)
            if chord_conductance:
                out = current.as_numpy()/(voltage-E_rev)
            else:
                out = abs(current.as_numpy())
            max_current[level] = max(out)
        return max_current
