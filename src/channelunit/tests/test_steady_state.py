import os
from collections import OrderedDict
import matplotlib.pyplot as plt
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
                 base_directory, save_figures):

        observation = self.format_data(observation)
        Test.__init__(self, observation, name)
        self.required_capabilities += (NModlChannel,)
        self.base_directory = base_directory
        self.save_figures = save_figures
        self.score_type = ZScore_SteadyStateCurves
        self.dpi = 200

    def compute_score(self, observation, prediction, verbose=False):
        score_avg, errors = ZScore_SteadyStateCurves.compute(observation,
                                                             prediction)
        score = ZScore_SteadyStateCurves(score_avg)
        return score

    def generate_figures(self, model, observations, predictions, name):
        v_values = list(observations.keys())
        pred_val = [predictions[v] for v in v_values]
        obs_val = [observations[v][0] for v in v_values]
        obs_std = [observations[v][1] for v in v_values]
        if self.base_directory:
            path = os.path.join(self.base_directory, "figs",
                                model.channel_name)
        else:
            path = os.path.join(model.base_directory, "figs",
                                model.channel_name)
        if not os.path.exists(path):
            os.makedirs(path)

        fig, ax = plt.subplots(1, 1)
        ax.plot(v_values, pred_val, "d", label=model.channel_name)
        ax.errorbar(v_values, obs_val, yerr=obs_std, marker="d", linewidth=0,
                    label="experimental data")
        ax.set_ylabel("Voltage (mV)")
        ax.set_xlabel("Normalized current")
        ax.legend()
        savefig_path = os.path.join(path, "%s_%s.png" % (model.channel_name,
                                                         name))
        fig.savefig(savefig_path, dpi=self.dpi,
                    bbox_inches='tight')

        
class ActivationSteadyStateTest(SteadyStateTest):

    def __init__(self, observation, experimental_conditions,
                 name="Activation Steady State Test",
                 base_directory="", save_figures=True):

        SteadyStateTest.__init__(self, observation, experimental_conditions,
                                name, base_directory, save_figures)
        try:
            self.v_init = experimental_conditions["v_init"]
        except KeyError:
            raise SystemExit("Please provide v_init in experimental_conditions")
        try:
            self.t_stop = experimental_conditions["t_stop"]
        except:
            raise SystemExit("Please provide stim length (t_stop (ms)) in experimental_conditions")
        try:
            self.chord_conductance = experimental_conditions["chord_conductance"]
        except:
            raise SystemExit("chord_conductance not specified in experimental conditions")    
        self.observation = OrderedDict(sorted(self.observation.items()))
        self.stimulus_list = self.extract_stimulation(self.observation)

    def run_model(self, model, stim_list, v_init, t_stop, chord_conductance):
        return model.get_activation_steady_state(stim_list,
                                                 v_init, t_stop,
                                                 chord_conductance)
    

    def generate_prediction(self, model, verbose=False):
        
        prediction = self.run_model(model, self.stimulus_list, self.v_init,
                                    self.t_stop, self.chord_conductance)
        if self.save_figures:
            name = self.name.replace(" ", "_")
            self.generate_figures(model, self.observation, prediction,
                                  name)
        return prediction


class InactivationSteadyStateTest(SteadyStateTest):

    def __init__(self, observation, experimental_conditions,
                 name="Inctivation Steady State Test",
                 base_directory="", save_figures=True):

        SteadyStateTest.__init__(self, observation, experimental_conditions,
                                name, base_directory, save_figures)
        try:
            self.v_test = experimental_conditions["v_test"]
        except KeyError:
            raise SystemExit("Please provide v_test in experimental_conditions")
        try:
            self.t_test = experimental_conditions["t_test"]
        except KeyError:
            raise SystemExit("Please provide stim length (t_test (ms)) in experimental_conditions")
        try:
            self.chord_conductance = experimental_conditions["chord_conductance"]
        except KeyError:
            raise SystemExit("chord_conductance not specified in experimental conditions")    
        self.observation = OrderedDict(sorted(self.observation.items()))
        self.stimulus_list = self.extract_stimulation(self.observation)

    def run_model(self, model, stim_list, v_test, t_test, chord_conductance):
        return model.get_inactivation_steady_state(stim_list,
                                                   v_test, t_test,
                                                   chord_conductance)
    

    def generate_prediction(self, model, verbose=False):
        
        prediction = self.run_model(model, self.stimulus_list, self.v_test,
                                    self.t_test, self.chord_conductance)
        if self.save_figures:
            name = self.name.replace(" ", "_")
            self.generate_figures(model, self.observation, prediction,
                                  name)

        return prediction

