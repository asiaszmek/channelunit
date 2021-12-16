import sciunit


class NModlChannel(sciunit.Capability):
    """This is an abstract model class providing capabilities
    that must be implemented by the patch model class.
    """
    def get_activation_steady_state(self, stimulation_levels: list,
                                    v_init: float, t_stop:float,
                                    chord_conductance=False,
                                    duration=200, sim_dt=0.001):
        """This function must be implemented by the patch model class.
        """

        raise NotImplementedError()

    def get_inactivation_steady_state(self, stimulation_levels: list,
                                      v_test: float, t_test:float,
                                      chord_conductance=False, sim_dt=0.001):
        """This function must be implemented by the patch model class.
        """

        raise NotImplementedError()

    def get_activation_traces(self, stimulation_levels: list,
                              v_hold: float, t_stop:float,
                              chord_conductance=False,
                              leak_subtraction=True,
                              duration=200, sim_dt=0.001, interval=200):
        """This function must be implemented by the patch model class.
        """

        raise NotImplementedError()
  
    def get_inactivation_traces(self, stimulation_levels: list,
                                v_test: float, t_test:float,
                                chord_conductance=False,
                                leak_subtraction=True,
                                sim_dt=0.001, interval=200):
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
        sim_dt: float
          for channels that can not be simulated using cvode. T
          his value should be small.
        """
        if self.cvode:
            h.cvode_active(1)
        else:
            h.cvode_active(0)
            h.dt = sim_dt

        h.celsius = self.temperature
        delay = 200
        stim_start = int(delay/self.dt)

        current_values = {}
        current = h.Vector()
        current.record(self.vclamp._ref_i, self.dt)
        time = h.Vector()
        time.record(h._ref_t, self.dt)
        for v_hold in stimulation_levels:
            t_stop = self.set_vclamp(delay, v_hold, t_test, v_test,
                                     leak_subtraction,
                                     delay=interval)
            h.init()
            h.tstop = t_stop
            h.run()
            I = current.as_numpy()[stim_start:]
            out = self.extract_current(I, chord_conductance,
                                       leak_subtraction, delay, t_test,
                                       interval)
            current_values[v_hold] = out
        return current_values
