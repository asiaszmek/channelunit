import os
import unittest

import numpy as np

from channelunit.tests import ActivationSteadyStateTest
from channelunit import ModelWholeCellPatchConcentration
from channelunit import data_path


channel_loc = os.path.join(data_path, "ion_channels")
activation_loc_Cal12 = os.path.join(data_path, "data",
                                    "I_CaL1.2_activation_Ca_1.5_mM.csv")
activation_loc_Cal13 = os.path.join(data_path, "data",
                                    "I_CaL1.3_activation_Ca_1.5_mM.csv")
activation_loc_Cal_20_Ba = os.path.join(data_path, "data",
                                        "I_CaL_activation_Ba_20_mM.csv")
activation_loc_Cal_110_Ba = os.path.join(data_path, "data",
                                         "I_CaL_activation_Ba_20_mM.csv")


class TestCaLChannelsLowBarium(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.modelBa20_H = ModelWholeCellPatchConcentration(channel_loc,
                                                           "calHGHK",
                                                           "ba",
                                                           external_conc=20,
                                                           temp=22,
                                                           liquid_junction_pot=0)
        cls.modelBa20_L = ModelWholeCellPatchConcentration(channel_loc,
                                                           "calGHK", "ba",
                                                           external_conc=20, temp=22,
                                                           liquid_junction_pot=0)
        cls.modelBa20_12 = ModelWholeCellPatchConcentration(channel_loc,
                                                            "cal12", "ba",
                                                            external_conc=20, temp=22,
                                                            liquid_junction_pot=0)
        cls.modelBa20_13 = ModelWholeCellPatchConcentration(channel_loc,
                                                            "cal13", "ba",
                                                            external_conc=20, temp=22,
                                                            liquid_junction_pot=0)

        activation_data = np.loadtxt(activation_loc_Cal_20_Ba, skiprows=1,
                                     delimiter=",")
        cls.power = 1
        cls.activation_data = dict(val.tolist() for val in activation_data)
        cls.test_Ba20 = ActivationSteadyStateTest(cls.activation_data,
                                                {"v_init": -90, "t_stop": 70,
                                                 "electrode_current": False,
                                                 "chord_conductance":False}, 1,
                                                "ActvationSSTest",
                                                save_figures=True)



    def test_summarize_H(self):
        self.score = self.test_Ba20.judge(self.modelBa20_H)
        self.score.summarize()

    def test_summarize_L(self):
        self.score = self.test_Ba20.judge(self.modelBa20_L)
        self.score.summarize()

    def test_summarize_12(self):
        self.score = self.test_Ba20.judge(self.modelBa20_12)
        self.score.summarize()

    def test_summarize_13(self):
        self.score = self.test_Ba20.judge(self.modelBa20_13)
        self.score.summarize()

    
    def test_run_model_H_keys(self):
        out = self.test_Ba20.run_model(self.modelBa20_H, self.test_Ba20.stimulus_list,
                                       self.test_Ba20.v_init,
                                       self.test_Ba20.t_stop,
                                       self.power,
                                       self.test_Ba20.chord_conductance,
                                       self.test_Ba20.electrode_current)
        self.assertEqual(list(out.keys()), self.test_Ba20.stimulus_list)

    def test_run_model_H_values(self):
        out = self.test_Ba20.run_model(self.modelBa20_H, self.test_Ba20.stimulus_list,
                                       self.test_Ba20.v_init,
                                       self.test_Ba20.t_stop,
                                       self.power,
                                       self.test_Ba20.chord_conductance,
                                       self.test_Ba20.electrode_current)
        values = np.array(list(out.values()))
        is_all_less_1 = np.all((values<=1))
        self.assertTrue(is_all_less_1)

    def test_run_model_L_keys(self):
        out = self.test_Ba20.run_model(self.modelBa20_L, self.test_Ba20.stimulus_list,
                                       self.test_Ba20.v_init, self.test_Ba20.t_stop,
                                       self.power,
                                       self.test_Ba20.chord_conductance,
                                       self.test_Ba20.electrode_current)
        self.assertEqual(list(out.keys()), self.test_Ba20.stimulus_list)

    def test_run_model_L_values(self):
        out = self.test_Ba20.run_model(self.modelBa20_L, self.test_Ba20.stimulus_list,
                                       self.test_Ba20.v_init,
                                       self.test_Ba20.t_stop,
                                       self.power,
                                       self.test_Ba20.chord_conductance,
                                       self.test_Ba20.electrode_current)
        values = np.array(list(out.values()))
        is_all_less_1 = np.all((values<=1))
        self.assertTrue(is_all_less_1)

        
    def test_run_model_12_keys(self):
        out = self.test_Ba20.run_model(self.modelBa20_12, self.test_Ba20.stimulus_list,
                                       self.test_Ba20.v_init,
                                       self.test_Ba20.t_stop,
                                       self.power,
                                       self.test_Ba20.chord_conductance,
                                       self.test_Ba20.electrode_current)
        self.assertEqual(list(out.keys()), self.test_Ba20.stimulus_list)

    def test_run_model_12_values(self):
        out = self.test_Ba20.run_model(self.modelBa20_12, self.test_Ba20.stimulus_list,
                                       self.test_Ba20.v_init,
                                       self.test_Ba20.t_stop,
                                       self.power,
                                       self.test_Ba20.chord_conductance,
                                       self.test_Ba20.electrode_current)
        values = np.array(list(out.values()))
        is_all_less_1 = np.all((values<=1))
        self.assertTrue(is_all_less_1)

    def test_run_model_13_keys(self):
        out = self.test_Ba20.run_model(self.modelBa20_13, self.test_Ba20.stimulus_list,
                                       self.test_Ba20.v_init, self.test_Ba20.t_stop,
                                       self.power,
                                       self.test_Ba20.chord_conductance,
                                       self.test_Ba20.electrode_current)
        self.assertEqual(list(out.keys()), self.test_Ba20.stimulus_list)

    def test_run_model_13_values(self):
        out = self.test_Ba20.run_model(self.modelBa20_13, self.test_Ba20.stimulus_list,
                                       self.test_Ba20.v_init,
                                       self.test_Ba20.t_stop,
                                       self.power,
                                       self.test_Ba20.chord_conductance,
                                       self.test_Ba20.electrode_current)
        values = np.array(list(out.values()))
        is_all_less_1 = np.all((values<=1))
        self.assertTrue(is_all_less_1)



if __name__ == "__main__":
    unittest.main()
