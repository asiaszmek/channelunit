import os
import unittest

import numpy as np

from channelunit.tests import SteadyStateTest
from channelunit import ModelWholeCellPatchSingleChan
from channelunit import data_path


channel_loc = os.path.join(data_path, "ion_channels")
activation_loc = os.path.join(data_path, "data",
                              "I_Na_activation.csv")
inactivation_loc = os.path.join(data_path, "data",
                                "I_Na_inactivation.csv")

class TestSteadyState(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.model = ModelWholeCellPatchSingleChan(channel_loc, "na3", "na",
                                                  external_conc=110,
                                                  temp=22, v_rest=-70,
                                                  ljp=0)
        
        activation_data = np.loadtxt(activation_loc, skiprows=1,
                                     delimiter=",")
        inactivation_data = np.loadtxt(inactivation_loc, skiprows=1,
                                       delimiter=",")
        cls.activation_data = dict(val.tolist() for val in activation_data)
         
        cls.inactivation_data = dict(val.tolist() for val in inactivation_data)
        
        cls.data = {"Activation": cls.activation_data,
                    "Inactivation": cls.inactivation_data}
        cls.experimental_conditions = {"Activation": {"v_init": -90,
                                                      "t_stop": 200,
                                                      "chord_conductance":
                                                      False},
                                       "Inactivation": {"v_test": -5,
                                                        "t_test": 10,
                                                        "chord_conductance":
                                                        False}
        }
        
        
        
        cls.test = SteadyStateTest(cls.data, cls.experimental_conditions,
                                   "Na3SS", electrode_current=False,
                                   normalization="to_one")
        
    def test_summarize(self):
        self.score = self.test.judge(self.model)
        self.score.summarize()


    def test_run_model(self):
        out = self.test.run_model(self.model, self.test.act_test.stimulus_list,
                                    self.test.act_test.v_init,
                                    self.test.act_test.t_stop,
                                    self.test.act_test.power,
                                    self.test.act_test.chord_conductance,
                                    self.test.inact_test.stimulus_list,
                                    self.test.inact_test.v_test,
                                    self.test.inact_test.t_test,
                                    self.test.inact_test.power,
                                    self.test.inact_test.chord_conductance,
                                    self.test.electrode_current,
                                    self.test.normalization)
        self.assertEqual(list(out.keys()), ["Activation", "Inactivation"])

        
if __name__ == "__main__":
    unittest.main()
