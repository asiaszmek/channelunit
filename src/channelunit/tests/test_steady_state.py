import sciunit
from sciunit import Test, Score

from channelunit.capabilities import NModlChannel


class SteadyStateTest(Test):
    """
    Common features of the Activation Steady State and Inactivation steady state
    """
    score_type = scores.ZScore_SteadyStateCurves
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
            elif isinstance(key, float):
                stim_list.append(key)
            elif isinstance(key, int):
                stim_list.append(key)
        return stim_list

    def format_data(self, observation_dict):
        pass
    
class ActivationSteadyStateTest(SteadyStateTest):

    def __init__(self, observation, name="Activation Steady State Test",
                 base_directory="", show_figures=True):

        self.injections = self.extract_stimulation(observation)
        observation = self.format_data(observation)
        observation = self.add_std_to_observation(observation)
        Test.__init__(self, observation, name)
        self.required_capabilities += (NModlChannel)
        self.base_directory = base_directory
        self.show_figures = show_figures
        self.npool = multiprocessing.cpu_count() - 1

    def format_data(self, observation):
        new_observation = {}
        for (key, val) in observation:
            if isinstance(str, val):
                values = val.split(" ")
                new_observation[key] = values[0]
                if len(values) == 2:
                    key_std = "std_" + str(key)
                    new_observation[key_std] = 
