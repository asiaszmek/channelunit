import os
from collections import OrderedDict
import matplotlib.pyplot as plt
from sciunit import Test, Score


from channelunit.capabilities import NModlChannel
from channelunit.scores import ZScore_SteadyStateCurves

class BaseSteadyStateTest(Test):
    """
    Common features of the Activation Steady State and Inactivation 
    steady state
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


    def add_to_fig(v_values_1, pred_val_1, obs_val_1, obs_std_1,
                   ax_1, label, marker="d",
                   xlabel="Voltage (mV)",
                   ylabel="Normalized current"):
        ax.plot(v_values_1, pred_val_1, marker, label=label)
        ax.errorbar(v_values_1, obs_val_1, yerr=obs_std_1,
                    marker=marker_1,
                    linewidth=0, label="experimental data")
        ax.set_xlabel(xlabel)
        ax.set_ylabel(ylabel)
        ax.legend()

    def generate_figures(self, model, observations, predictions, name):
        v_values = list(observations.keys())
        pred_val = [predictions[v] for v in v_values]
        obs_val = [observations[v][0] for v in v_values]
        obs_std = [observations[v][1] for v in v_values]
        channel_names_path = ""
        for channel_name in model.channel_names:
            channel_names_path += channel_name + "_"
        if self.base_directory:
            path = os.path.join(self.base_directory, "figs",
                                channel_names_path)
        else:
            path = os.path.join(model.base_directory, "figs",
                                channel_names_path)
        if not os.path.exists(path):
            os.makedirs(path)

        fig, ax = plt.subplots(1, 1)
        label = ""
        for n in model.channel_names:
            label += n + " "
        self.add_to_fig(v_values, pred_val, obs_val, obs_std,
                        ax, label)
        savefig_path = os.path.join(path, "%s_%s.png" % (channel_names_path,
                                                         name))
        fig.savefig(savefig_path, dpi=self.dpi,
                    bbox_inches='tight')

        
class ActivationSteadyStateTest(BaseSteadyStateTest):

    def __init__(self, observation, experimental_conditions, power: int,
                 name="Activation Steady State Test",
                 base_directory="", save_figures=True):
        conditions = experimental_conditions
        super(ActivationSteadyStateTest, self).__init__(observation,
                                                            conditions,
                                                            name,
                                                            base_directory,
                                                            save_figures)
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
        try:
            self.electrode_current = experimental_conditions["electrode_current"]
        except KeyError:
            raise SystemExit("It must be specified, if it's an electrode_current")
        try:
            self.normalization = experimental_conditions["normalization"]
        except KeyError:
            raise SystemExit("It must be specified, if currents normalized to_one")    
        self.power = power
        self.observation = OrderedDict(sorted(self.observation.items()))
        self.stimulus_list = self.extract_stimulation(self.observation)

    def run_model(self, model, stim_list, v_init, t_stop,
                  power, chord_conductance, electrode_current, normalization):
        return model.get_activation_SS(stim_list,
                                       v_init, t_stop, power,
                                       chord_conductance,
                                       electrode_current,
                                       normalization=normalization)
    

    def generate_prediction(self, model, verbose=False):
        
        prediction = self.run_model(model, self.stimulus_list, self.v_init,
                                    self.t_stop, self.power,
                                    self.chord_conductance,
                                    self.electrode_current,
                                    self.normalization)
        if self.save_figures:
            name = self.name.replace(" ", "_")
            self.generate_figures(model, self.observation, prediction,
                                  name)
        return prediction


class InactivationSteadyStateTest(BaseSteadyStateTest):

    def __init__(self, observation, experimental_conditions, power,
                 name="Inctivation Steady State Test",
                 base_directory="", save_figures=True):
        conditions = experimental_conditions
        super(InactivationSteadyStateTest, self).__init__(observation,
                                                              conditions,
                                                              name,
                                                              base_directory,
                                                              save_figures)
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
        try:
            self.electrode_current = experimental_conditions["electrode_current"]
        except KeyError:
              raise SystemExit("It must be specified, if it's an electrode_current")
        try:
            self.normalization = experimental_conditions["normalization"]
        except KeyError:
            raise SystemExit("It must be specified, if currents normalized to_one")    
        self.power = power
        self.observation = OrderedDict(sorted(self.observation.items()))
        self.stimulus_list = self.extract_stimulation(self.observation)

    def run_model(self, model, stim_list, v_test, t_test, power,
                  chord_conductance, electrode_current, normalization):
        return model.get_inactivation_SS(stim_list,
                                         v_test, t_test, power,
                                         chord_conductance,
                                         electrode_current,
                                         normalization=normalization)
    

    def generate_prediction(self, model, verbose=False):
        
        prediction = self.run_model(model, self.stimulus_list, self.v_test,
                                    self.t_test, self.power,
                                    self.chord_conductance,
                                    self.electrode_current,
                                    self.normalization)
        if self.save_figures:
            name = self.name.replace(" ", "_")
            self.generate_figures(model, self.observation, prediction,
                                  name)

        return prediction

