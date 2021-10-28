import sciunit
from sciunit import Test, Score

from channelunit.capabilities import NModlChannel


class SteadyStateTest(Test):
    """
    Common features of the Activation Steady State and Inactivation steady state
    """
    #score_type = scores.ZScore_SteadyStateCurves
    def extract_stimulation(self, observation_dict):
        stim_list = []
        for key in observation_dict.keys():
            if isinstance(key, str):
                try:
                    numerical_key = float(key)
                except ValueError:
                    print("observation %s is not a float" % key)
                    continue
                stim_list.append(float(key))
            elif isinstance(key, float) or isinstance(key, int):
                stim_list.append(key)
        return stim_list

    def format_data(self, observation_dict):
        observation_std = {}
        observation = {}
        for (key, val) in observation_dict:
            if isinstance(key, str):
                new_key = float(key)
            elif isinstance(key, float):
                new_key = key
            elif isinstance(key, int):
                new_key = key
            else:
                continue
            if isinstance(val, str):
                values = val.split(" ")
                new_val = float(values[0])
                if len(values) == 2:
                    std_val = float(values[1])
                    observation_std[new_key] = std_val
            elif isinstance(val, float) or isinstance(val, int):
                new_val = val
            observation[new_key] = new_val
            if new_key not in observation_std:
                observation_std[new_key] = 0.03
        return observation, observation_std

    
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

       
