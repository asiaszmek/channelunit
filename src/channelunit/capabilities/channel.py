import sciunit


class NModlChannel(sciunit.Capability):
    """This is an abstract model class providing capabilities
    that must be implemented by the patch model class.
    """
    def get_activation_SS(self, stimulation_levels: list,
                          v_hold: float, t_stop:float,
                          power: int, t_mes, chord_conductance,
                          electrode_current,
                          normalization="to_one",
                          save_traces=True, save_ca=False):
        """This function must be implemented by the patch model class.
        """

        raise NotImplementedError()

    def get_inactivation_SS(self, stimulation_levels: list,
                            v_test: float, t_test:float,
                            power: int, t_mes,
                            chord_conductance=False,
                            electrode_current=True,
                            normalization="to_one",
                            save_traces=True, save_ca=True):
        """This function must be implemented by the patch model class.
        """

        raise NotImplementedError()

    def get_activation_traces(self, stimulation_levels: list,
                              v_hold: float, t_stop:float,
                              chord_conductance=False,
                              electrode_current=True,
                              save_traces=True, save_ca=False):
        """This function must be implemented by the patch model class.
        """

        raise NotImplementedError()
  
    def get_inactivation_traces(self, stimulation_levels: list,
                                v_test: float, t_test:float,
                                chord_conductance,
                                electrode_current,
                                save_traces=True, save_ca=True):

        
        """This function must be implemented by the patch model class.
        """

        raise NotImplementedError()
