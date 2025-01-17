import os
import unittest

import numpy as np

from channelunit.tests import ActivationSteadyStateTest
from channelunit import ModelWholeCellPatchCaSingleChan
from channelunit import data_path


channel_loc = os.path.join(data_path, "ion_channels")

activation_loc_CaT_10_ca = os.path.join(data_path, "data",
                                           "TC_CaT_Stuart_act.csv")
inactivation_loc_CaT_10_ca = os.path.join(data_path, "data",
                                             "TC_CaT_Stuart_inact.csv")

class TestCaTChannels_ca(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.modelca_10_H = ModelWholeCellPatchCaSingleChan(channel_loc,
                                                          "caT",
                                                          "ca",
                                                           external_conc=10,
                                                           temp=25,
                                                           ljp=0,
                                                           gbar_value=0.001)
        activation_data = np.loadtxt(activation_loc_CaT_10_ca, skiprows=1,
                                     delimiter=",")
        cls.power = 2
        cls.activation_data = dict(val.tolist() for val in activation_data)
        cls.test_ca10 = ActivationSteadyStateTest(cls.activation_data,
                                                  {"v_init": -20, "t_stop": 120,
                                                   "electrode_current": False,
                                                   "chord_conductance":False,
                                                   "normalization": "to_one"},
                                                  1,
                                                  "ActvationSSTestCaTca",
                                                  save_figures=True)
        cls.act_results = cls.test_ca10.run_model(cls.modelca_10_H,
                                                  cls.test_ca10.stimulus_list,
                                                  cls.test_ca10.v_init,
                                                  cls.test_ca10.t_stop,
                                                  cls.power,
                                                  cls.test_ca10.chord_conductance,
                                                  cls.test_ca10.electrode_current,
                                                  "to_one")
       

    def test_summarize_H(self):
        self.score = self.test_ca10.judge(self.modelca_10_H)
        self.score.summarize()

    def test_run_model_H_keys(self):
        self.assertEqual(list(self.act_results.keys()),
                         self.test_ca10.stimulus_list)

    def test_run_model_H_values(self):
        values = np.array(list(self.act_results.values()))
        is_all_less_1 = np.all((values<=1))
        self.assertTrue(is_all_less_1)

    def test_gbar_val(self):
        gbar_val = self.modelca_10_H.patch.psection()["density_mechs"]["caT"]["gbar"]
        print(self.modelca_10_H.ca_decay)
        print(self.modelca_10_H.decay_eq)
        print(self.modelca_10_H.patch.psection())
        self.assertEqual(gbar_val, [0.0002])
        


if __name__ == "__main__":
    unittest.main()
