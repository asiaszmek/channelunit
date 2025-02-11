import os
import unittest

import numpy as np

from channelunit.tests import ActivationSteadyStateTest
from channelunit.tests import InactivationSteadyStateTest

from channelunit import ModelWholeCellPatchCaSingleChan
from channelunit import data_path


channel_loc = os.path.join(data_path, "ion_channels")

activation_loc_CaT_2_ca = os.path.join(data_path, "data",
                                       "McRory_cat_a1g_act.csv")
inactivation_loc_CaT_2_ca = os.path.join(data_path, "data",
                                         "McRory_cat_a1g_inact.csv")

class TestCaTChannels_ca_Act(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.modelca_2_H = ModelWholeCellPatchCaSingleChan(channel_loc,
                                                          "ca12dZUy",
                                                          "ca",
                                                           external_conc=2,
                                                           temp=22,
                                                           ljp=0, v_rest=-65,
                                                           gbar_value=0.001)
        activation_data = np.loadtxt(activation_loc_CaT_2_ca, skiprows=1,
                                     delimiter=",")
        cls.power = 1
        cls.activation_data = dict(val.tolist() for val in activation_data)
        cls.act_test_ca2 = ActivationSteadyStateTest(cls.activation_data,
                                                  {"v_init": -80, "t_stop": 15,
                                                   "electrode_current": False,
                                                   "chord_conductance":False,
                                                   "normalization": "to_one"},
                                                  1,
                                                  "ActivationSSTestCaTca",
                                                  save_figures=True)

        cls.act_results = cls.act_test_ca2.run_model(cls.modelca_2_H,
                                                     cls.act_test_ca2.stimulus_list,
                                                     cls.act_test_ca2.v_init,
                                                     cls.act_test_ca2.t_stop,
                                                     cls.power,
                                                     cls.act_test_ca2.t_mes,
                                                     cls.act_test_ca2.chord_conductance,
                                                     cls.act_test_ca2.electrode_current,
                                                  "to_one")

    def test_act_summarize_H(self):
        self.score = self.act_test_ca2.judge(self.modelca_2_H)
        self.score.summarize()

    def test_act_run_model_H_keys(self):
        self.assertEqual(list(self.act_results.keys()),
                         self.act_test_ca2.stimulus_list)

    def test_act_run_model_H_values(self):
        values = np.array(list(self.act_results.values()))
        is_all_less_1 = np.all((values<=1))
        self.assertTrue(is_all_less_1)



class TestCaTChannels_ca_Inact(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.modelca_2_H = ModelWholeCellPatchCaSingleChan(channel_loc,
                                                          "ca12dZUy",
                                                          "ca",
                                                           external_conc=2,
                                                           temp=22,
                                                           ljp=0, v_rest=-65,
                                                           gbar_value=0.001)
        inactivation_data = np.loadtxt(inactivation_loc_CaT_2_ca, skiprows=1,
                                        delimiter=",")
        cls.power = 1
        cls.inactivation_data = dict(val.tolist() for val in inactivation_data)
        
        cls.inact_test_ca2 = InactivationSteadyStateTest(cls.inactivation_data,
                                                  {"v_test": -30, "t_test": 1000,
                                                   "electrode_current": False,
                                                   "chord_conductance":False,
                                                   "normalization": "to_one"},
                                                  1,
                                                  "InactivationSSTestCaTca",
                                                  save_figures=True)
        cls.inact_results = cls.inact_test_ca2.run_model(cls.modelca_2_H,
                                                         cls.inact_test_ca2.stimulus_list,
                                                         cls.inact_test_ca2.v_test,
                                                         cls.inact_test_ca2.t_test,
                                                         cls.power,
                                                         cls.inact_test_ca2.t_mes,
                                                         cls.inact_test_ca2.chord_conductance,
                                                         cls.inact_test_ca2.electrode_current,
                                                         "to_one")
       

    def test_inact_summarize_H(self):
        self.score = self.inact_test_ca2.judge(self.modelca_2_H)
        self.score.summarize()

    def test_inact_run_model_H_keys(self):
        self.assertEqual(list(self.inact_results.keys()),
                         self.inact_test_ca2.stimulus_list)

    def test_inact_run_model_H_values(self):
        values = np.array(list(self.inact_results.values()))
        is_all_less_1 = np.all((values<=1))
        self.assertTrue(is_all_less_1)



if __name__ == "__main__":
    unittest.main()
