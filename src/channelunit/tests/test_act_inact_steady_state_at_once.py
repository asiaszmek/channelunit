import os
from collections import OrderedDict
import matplotlib.pyplot as plt
from sciunit import Test, Score


from channelunit.capabilities import NModlChannel
from .test_steady_state import InactivationSteadyStateTest
from .test_steady_state import BaseSteadyStateTest
from .test_steady_state import ActivationSteadyStateTest
from channelunit.scores import ZScore_BothSteadyStateCurves

class SteadyStateTest(Test):
    """
    Test if steady state activation curve and steady state inactivation
    curves of the model channel fit experimental data

    observation: a dictonary with "Activation" and "Inactivation" keys. 
          Values are dictionaries of experimental observations in the form:
          key: applied voltage, value: activation value
    """
    def __init__(self, observation: dict, experimental_conditions: dict,
                 name, electrode_current, normalization,
                 power={"Activation": 1, "Inactivation": 1},
                 base_directory="", save_figures=True):
        act_obs = observation["Activation"]
        act_cond = experimental_conditions["Activation"]
        act_cond["electrode_current"] = electrode_current
        act_cond["normalization"] = normalization
        inact_obs = observation["Inactivation"]
        inact_cond = experimental_conditions["Inactivation"]
        inact_cond["electrode_current"] = electrode_current
        inact_cond["normalization"] = normalization
        self.act_test = ActivationSteadyStateTest(act_obs, act_cond,
                                                  power["Activation"],
                                                  "%s_%s" %(name,
                                                            "Activation"),
                                                  base_directory,
                                                  save_figures=False)
        
        self.inact_test = InactivationSteadyStateTest(inact_obs, inact_cond,
                                                      power["Inactivation"],
                                                      "%s_%s" % (name,
                                                                 "Inactivation"),
                                                    base_directory,
                                                    save_figures=False)
        observation = {"Activation": self.act_test.observation,
                       "Inactivation": self.inact_test.observation}
        super(SteadyStateTest, self).__init__(observation, name)
        self.electrode_current = electrode_current
        self.power = power
        self.normalization = normalization
        self.required_capabilities += (NModlChannel,)
        self.base_directory = base_directory
        self.save_figures = save_figures
        self.score_type = ZScore_BothSteadyStateCurves
        self.dpi = 200
    
    def compute_score(self, observation, prediction, verbose=False):
        score_avg, errors = ZScore_BothSteadyStateCurves.compute(observation,
                                                                 prediction)
        score = ZScore_BothSteadyStateCurves(score_avg)
        return score

        
    def run_model(self, model, act_stim_list, act_v_init, act_t_stop,
                  act_power, act_chord_conductance, inact_stim_list,
                  inact_v_test, inact_t_test,
                  inact_power,  inact_chord_conductance,
                  electrode_current, normalization):
        act_prediction = model.get_activation_SS(act_stim_list,
                                                 act_v_init, act_t_stop,
                                                 act_power,
                                                 act_chord_conductance,
                                                 electrode_current,
                                                 normalization,
                                                 save_traces=False,
                                                 save_ca=False)
        inact_prediction = model.get_inactivation_SS(inact_stim_list,
                                                     inact_v_test,
                                                     inact_t_test,
                                                     inact_power,
                                                     inact_chord_conductance,
                                                     electrode_current,
                                                     normalization,
                                                     save_traces=False,
                                                     save_ca=False)
        return {"Activation": act_prediction,
                "Inactivation": inact_prediction}
    

    def generate_prediction(self, model, verbose=False):
        
        prediction = self.run_model(model, self.act_test.stimulus_list,
                                    self.act_test.v_init,
                                    self.act_test.t_stop,
                                    self.act_test.power,
                                    self.act_test.chord_conductance,
                                    self.inact_test.stimulus_list,
                                    self.inact_test.v_test,
                                    self.inact_test.t_test,
                                    self.inact_test.power,
                                    self.inact_test.chord_conductance,
                                    self.electrode_current,
                                    self.normalization)
        if self.save_figures:
            self.generate_figures(model, self.observation, prediction, "SST")
        return prediction


    def generate_figures(self, model, observations, predictions, name):
        label = ""
        for n in model.channel_names:
            label += n + " "
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
        
        act_v = list(observations["Activation"].keys())
        act_p = [predictions["Activation"][v] for v in act_v]
        act_vals = [observations["Activation"][v][0] for v in act_v]
        act_vals_std = [observations["Activation"][v][1] for v in act_v]
        if self.act_test.chord_conductance:
            label1 = "Normalized conductance"
        else:
            label1 = "Normalized current"
        if self.inact_test.chord_conductance:
            label2 = "Normalized conductance"
        else:
            label2 = "Normalized current"

        ax2 = ax.twinx()
        BaseSteadyStateTest.add_to_figure(act_v, act_p, act_vals, act_vals_std,
                                          ax2, label,
                                          ylabel=label1)
        
        
        inact_v = list(observations["Inactivation"].keys())
        inact_p = [predictions["Inactivation"][v] for v in inact_v]
        inact_vals = [observations["Inactivation"][v][0] for v in inact_v]
        inact_vals_std = [observations["Inactivation"][v][1] for v in inact_v]
        
        BaseSteadyStateTest.add_to_figure(inact_v, inact_p, inact_vals,
                                          inact_vals_std,
                                          ax, label, marker="o", color="k",
                                          ylabel=label2)
        ymax = max(max(ax.get_ylim()), max(ax2.get_ylim()))
        ymin = min(min(ax.get_ylim()), min(ax2.get_ylim()))
        ax.set_ylim([ymin, ymax])
        ax2.set_ylim([ymin, ymax])
        ax.legend(loc=3)
        savefig_path = os.path.join(path, "Steady_state_%s_%s.png"
                                    % (channel_names_path, name))
        fig.savefig(savefig_path, dpi=self.dpi,
                    bbox_inches='tight')
