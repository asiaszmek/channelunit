import sciunit
from sciunit import Test, Score

from channelunit.capabilities import NModlChannel
from channelunit.scores import ZScore_SteadyStateCurves

class SteadyStateTest(Test):
    """
    Common features of the Activation Steady State and Inactivation steady state
    """

    @classmethod
    def _get_numerical(self, val):
        if isinstance(val, str):
            try:
                numerical_key = float(val)
            except ValueError:
                raise ValueError("observation is not a float")
            return numerical_key
        if isinstance(val, int) or isinstance(val, float):
            return val
        if isinstance(val, list):
            if len(val) == 1:
                if isinstance(val[0], int) or isinstance(val[0], float):
                    return val[0]
                if isinstance(val[0], str):
                    try:
                        numerical_key = float(val[0])
                    except ValueError:
                        raise ValueError("observation is not a float")
                    return numerical_key

        raise ValueError("observation is not a float")

    @classmethod
    def extract_stimulation(self, observation_dict):
        stim_list = []
        for key in observation_dict.keys():
            stim_list.append(self._get_numerical(key))
        return sorted(stim_list)

    @classmethod
    def format_data(self, observation_dict):
        observation = {}
        for key in observation_dict.keys():
            new_key = self._get_numerical(key)
            val = observation_dict[key]
            if not isinstance(val, list):
                val = [val]
            new_val = []
            for v in val:
                new_val.append(self._get_numerical(v))
            if len(new_val) == 1:
                new_val.append(0.03)
            observation[new_key] = new_val
        return observation

class ActivationSteadyStateTest(SteadyStateTest):

    def __init__(self, observation, name="Activation Steady State Test",
                 base_directory="", show_figures=True):

        self.injections = self.extract_stimulation(observation) 
        observation_vals, observation_stds = self.format_data(observation)
        Test.__init__(self, observation, name)
        self.required_capabilities += (NModlChannel)
        self.base_directory = base_directory
        self.show_figures = show_figures
        self.npool = multiprocessing.cpu_count() - 1
        self.score_type = ZScore_SteadyStateCurves
