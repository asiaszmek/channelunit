from collections import OrderedDict
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

    def __init__(self, observation, experimental_conditions, name,
                 base_directory, show_figures):

        self.observation = self.format_data(observation)
        Test.__init__(self, observation, name)
        self.required_capabilities += (NModlChannel)
        self.base_directory = base_directory
        self.show_figures = show_figures
        self.score_type = ZScore_SteadyStateCurves
        self.experimental_conditions = experimental_conditions

    def compute_score(self, observation, prediction, verbose=False):
        score_avg, errors = ZScore_SteadyStateCurves(observation,
                                                     prediction)
        return score_avg
        

class ActivationSteadyStateTest(SteadyStateTest):

    def __init__(self, observation, experimental_conditions,
                 name="Activation Steady State Test",
                 base_directory="", show_figures=True):

        SteadyStateTest._init__(self, observation, experimental_conditions,
                                name, base_directory, show_figures)
        if "v_init" not in self.experimental_conditions:
            raise SystemExit("Please provide v_init in experimental_conditions")
        if "t_stop" not in self.experimental_conditions:
            raise SystemExit("Please provide stim length (t_stop (ms)) in experimental_conditions")
        if chord_conductance not in self.experimental_conditions:
            raise SystemExit("chord_conductance not specified in experimental conditions")    

    def run_model(self, model, stim_list, v_init, t_stop, chord_conductance):
        return model.get_activation_steady_state(stim_list,
                                                 v_init, t_stop,
                                                 chord_conductance)
    

    def generate_prediction(self, model, verbose=False):
        
        self.observation = OrderedDict(sorted(self.observation.items()))
        self.stimulus_list = self.extract_stimulation(self.observation)
        v_init = self.experimental_conditions["v_init"]
        t_stop = self.experimental_conditions["t_stop"]
        chord_conductance = self.expermental_conditions["chord_conductance"]
        prediction = self.run_model(model, self.stimulus_list, v_init, t_stop,
                                    chord_conductance)
        return prediction


class InactivationSteadyStateTest(SteadyStateTest):

    def __init__(self, observation, experimental_conditions,
                 name="Inctivation Steady State Test",
                 base_directory="", show_figures=True):

        SteadyStateTest._init__(self, observation, experimental_conditions,
                                name, base_directory, show_figures)
        if "v_test" not in self.experimental_conditions:
            raise SystemExit("Please provide v_test in experimental_conditions")
        if "t_test" not in self.experimental_conditions:
            raise SystemExit("Please provide stim length (t_test (ms)) in experimental_conditions")
        if chord_conductance not in self.experimental_conditions:
            raise SystemExit("chord_conductance not specified in experimental conditions")    

    def run_model(self, model, stim_list, v_test, t_test, chord_conductance):
        return model.get_inactivation_steady_state(stim_list,
                                                   v_test, t_test,
                                                   chord_conductance)
    

    def generate_prediction(self, model, verbose=False):
        
        self.observation = OrderedDict(sorted(self.observation.items()))
        self.stimulus_list = self.extract_stimulation(self.observation)
        v_test = self.experimental_conditions["v_test"]
        t_test = self.experimental_conditions["t_test"]
        chord_conductance = self.expermental_conditions["chord_conductance"]
        prediction = self.run_model(model, self.stimulus_list, v_test,
                                    t_test, chord_conductance)
        return prediction

