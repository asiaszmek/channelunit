import os
import unittest

import numpy as np

from channelunit.tests import ActivationSteadyStateTest
from channelunit import ModelWholeCellPatchSingleChan
from channelunit import data_path


channel_loc = os.path.join(data_path, "ion_channels")
# proximal A-type K channels
activation_loc_K_A_p = os.path.join(data_path, "data",
                                    "I_K_A_prox_activation.csv")

#distal A-type K channels
activation_loc_K_A_d = os.path.join(data_path, "data",
                                    "I_K_A_dist_activation.csv")
                                   
inactivation_loc_K_A = os.path.join(data_path, "data",
                                    "I_K_A_inactivation.csv")
                       

class TestProximalKA(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.modelKad = ModelWholeCellPatchSingleChan(channel_loc,
                                                     "kad",
                                                     "k",
                                                     external_conc=2.5,
                                                     temp=22,
                                                     ljp=0,
                                                     v_rest=-65)
        activation_data = np.loadtxt(activation_loc_K_A_d, skiprows=1,
                                     delimiter=",")
        inactivation_data = np.loadtxt(inactivation_loc_K_A, skiprows=1,
                                       delimiter=",")

        cls.power = 1
        cls.activation_data = dict()
        for val in activation_data:
            cls.activation_data[val[0]] = val[1:].tolist()
        cls.inactivation_data = dict()
        for val in inactivation_data:
            cls.inactivation_data[val[0]] = val[1:].tolist()
            
        cls.test_K_A_type = ActivationSteadyStateTest(cls.activation_data,
                                                      {"v_init": -85, "t_stop": 50,
                                                       "electrode_current": True,
                                                       "chord_conductance": True,
                                                       "normalization":
                                                       "save_sign"}, 1,
                                                      "ActvationSSTest",
                                                      save_figures=True)
        
        cls.act_results = cls.test_K_A_type.run_model(cls.modelKad,
                                                      cls.test_K_A_type.stimulus_list,
                                                      cls.test_K_A_type.v_init,
                                                      cls.test_K_A_type.t_stop,
                                                      cls.power,
                                                      cls.test_K_A_type.chord_conductance,
                                                      cls.test_K_A_type.electrode_current,
                                                      "save_sign")
       

    def test_summarize(self):
        self.score = self.test_K_A_type.judge(self.modelKad)
        self.score.summarize()

    def test_run_model_keys(self):
        self.assertEqual(list(self.act_results.keys()),
                         self.test_K_A_type.stimulus_list)

    def test_run_model_values(self):
        values = np.array(list(self.act_results.values()))
        is_all_less_1 = np.all((values<=1))
        self.assertTrue(is_all_less_1)

        
class TestDistalKA(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        
        cls.modelKap = ModelWholeCellPatchSingleChan(channel_loc,
                                                     "kap",
                                                     "k",
                                                     external_conc=2.5,
                                                     temp=22,
                                                     ljp=0)

        activation_data = np.loadtxt(activation_loc_K_A_p, skiprows=1,
                                     delimiter=",")
        inactivation_data = np.loadtxt(inactivation_loc_K_A, skiprows=1,
                                       delimiter=",")

        cls.power = 1
        cls.activation_data = dict()
        for val in activation_data:
            cls.activation_data[val[0]] = val[1:].tolist()
        cls.test_K_A_type = ActivationSteadyStateTest(cls.activation_data,
                                                {"v_init": -85, "t_stop": 50,
                                                 "electrode_current": False,
                                                 "chord_conductance": True,
                                                 "normalization": "save_sign"}, 1,
                                                "ActvationSSTest",
                                                save_figures=True)
        cls.act_results = cls.test_K_A_type.run_model(cls.modelKap,
                                                  cls.test_K_A_type.stimulus_list,
                                                  cls.test_K_A_type.v_init,
                                                  cls.test_K_A_type.t_stop,
                                                  cls.power,
                                                  cls.test_K_A_type.chord_conductance,
                                                  cls.test_K_A_type.electrode_current,
                                                  "save_sign")


    def test_summarize(self):
        self.score = self.test_K_A_type.judge(self.modelKap)
        self.score.summarize()

    def test_run_model_keys(self):
        self.assertEqual(list(self.act_results.keys()),
                         self.test_K_A_type.stimulus_list)

    def test_run_model_values(self):
        values = np.array(list(self.act_results.values()))
        
        is_all_less_1 = np.all((values<=1))

        self.assertTrue(is_all_less_1)


if __name__ == "__main__":
    unittest.main()
