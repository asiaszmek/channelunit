import sciunit


class NModlChannel(sciunit.Capability):
    """This is an abstract model class providing capabilities
    that must be implemented by the patch model class.
    """
    def get_activation_SS(self, stimulation_levels: list,
                          v_init: float, t_stop:float,
                          chord_conductance=False,
                          electrode_current=True,
                          sim_dt=0.001, normalization="to_one"):
        """This function must be implemented by the patch model class.
        """

        raise NotImplementedError()

    def get_inactivation_SS(self, stimulation_levels: list,
                            v_test: float, t_test:float,
                            chord_conductance=False,
                            electrode_current=True,
                            sim_dt=0.001, normalization="to_one"):
        """This function must be implemented by the patch model class.
        """

        raise NotImplementedError()

    def get_activation_traces(self, stimulation_levels: list,
                              v_hold: float, t_stop:float,
                              chord_conductance=False,
                              electrode_current=True,
                              sim_dt=0.001, interval=200):
        """This function must be implemented by the patch model class.
        """

        raise NotImplementedError()
  
    def get_inactivation_traces(self, stimulation_levels: list,
                                v_test: float, t_test:float,
                                chord_conductance=False,
                                electrode_current=True,
                                sim_dt=0.001, interval=200):
        
        """This function must be implemented by the patch model class.
        """

        raise NotImplementedError()
