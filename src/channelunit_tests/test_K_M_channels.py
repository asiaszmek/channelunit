import os
import unittest

import numpy as np

from channelunit.tests import ActivationSteadyStateTest
from channelunit import ModelWholeCellPatchSingleChan
from channelunit import data_path


channel_loc = os.path.join(data_path, "ion_channels")
# proximal A-type K channels
activation_loc_K_M = os.path.join(data_path, "data",
                                  "IM_neocortical_axons.csv")

                       

class TestProximalKM(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.modelKM = ModelWholeCellPatchSingleChan(channel_loc,
                                                    "iM",
                                                    "k",
                                                    external_conc=2.5,
                                                    temp=35,
                                                    ljp=0,
                                                    v_rest=-106)

        activation_data = np.loadtxt(activation_loc_K_M, skiprows=1,
                                     delimiter=",")
        

        cls.power = 1
        cls.activation_data = dict()
        for val in activation_data:
            cls.activation_data[val[0]] = val[1:].tolist()
            
        cls.test_K_M_type = ActivationSteadyStateTest(cls.activation_data,
                                                      {"v_init": -85, "t_stop": 100,
                                                       "electrode_current": True,
                                                       "chord_conductance": True,
                                                       "normalization":
                                                       "save_sign"}, 1,
                                                      "ActvationSSTest",
                                                      save_figures=True)
        
        cls.act_results = cls.test_K_M_type.run_model(cls.modelKM,
                                                      cls.test_K_M_type.stimulus_list,
                                                      cls.test_K_M_type.v_init,
                                                      cls.test_K_M_type.t_stop,
                                                      cls.power,
                                                      cls.test_K_M_type.t_mes,
                                                      cls.test_K_M_type.chord_conductance,
                                                      cls.test_K_M_type.electrode_current,
                                                      "save_sign")
       

    def test_summarize(self):
        self.score = self.test_K_M_type.judge(self.modelKM)
        self.score.summarize()

    def test_run_model_keys(self):
        self.assertEqual(list(self.act_results.keys()),
                         self.test_K_M_type.stimulus_list)

    def test_run_model_values(self):
        values = np.array(list(self.act_results.values()))
        is_all_less_1 = np.all((values<=1))
        self.assertTrue(is_all_less_1)

        

if __name__ == "__main__":
    unittest.main()
