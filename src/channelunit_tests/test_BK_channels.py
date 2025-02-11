import os
import unittest

import numpy as np

from channelunit.tests import ActivationSteadyStateTest
from channelunit import ModelCaConcClamp
from channelunit import data_path


channel_loc = os.path.join(data_path, "ion_channels")
activation_loc_BK_1_uM_Ca = os.path.join(data_path, "data",
                                         "I_BKCa_1_uM_Ca.csv")
activation_loc_BK_10_uM_Ca = os.path.join(data_path, "data",
                                         "I_BKCa_10_uM_Ca.csv")
activation_loc_BK_100_uM_Ca = os.path.join(data_path, "data",
                                         "I_BKCa_100_uM_Ca.csv")

class TestBK_1_uM_Ca(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.model = ModelCaConcClamp(channel_loc,
                                     ["bk"],
                                     ["k", "ca"],
                                     gbar_values={"bk": 1},
                                     external_conc={"k": 5},
                                     internal_conc={"ca": 1e-3,
                                                    "k": 120},
                                     v_rest=-100,
                                     temp=23,
                                     ljp=0)
        activation_data = np.loadtxt(activation_loc_BK_1_uM_Ca, skiprows=1,
                                     delimiter=",")
        cls.activation_data = {}
        for val in activation_data:
            cls.activation_data[val[0]] = val[1:].tolist()
        cls.power = 1
        cls.test = ActivationSteadyStateTest(cls.activation_data,
                                             {"v_init": -80, "t_stop": 40,
                                              "electrode_current": True,
                                              "chord_conductance":False,
                                              "normalization": "to_one"},
                                             1,
                                             "ActvationTestBKCa1uMCa",
                                             save_figures=True)
        cls.act_results = cls.test.run_model(cls.model,
                                             cls.test.stimulus_list,
                                             cls.test.v_init,
                                             cls.test.t_stop,
                                             cls.power,
                                             cls.test.t_mes,
                                             cls.test.chord_conductance,
                                             cls.test.electrode_current,
                                             "to_one")

    def test_summarize(self):
        self.score = self.test.judge(self.model)
        self.score.summarize()

    def test_run_model_keys(self):
        self.assertEqual(list(self.act_results.keys()),
                         self.test.stimulus_list)

    def test_run_model_values(self):
        values = np.array(list(self.act_results.values()))
        is_all_less_1 = np.all((values<=1))
        self.assertTrue(is_all_less_1)


class TestBK_10_uM_Ca(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.model = ModelCaConcClamp(channel_loc,
                                     ["bk"],
                                     ["k", "ca"],
                                     gbar_values={"bk": 1},
                                     external_conc={"k": 5},
                                     internal_conc={"ca": 10e-3,
                                                    "k": 120},
                                     temp=23,
                                     v_rest=-100,
                                     ljp=0)
        activation_data = np.loadtxt(activation_loc_BK_10_uM_Ca, skiprows=1,
                                     delimiter=",")
        cls.activation_data = {}
        for val in activation_data:
            cls.activation_data[val[0]] = val[1:].tolist()
        cls.power = 1
        cls.test = ActivationSteadyStateTest(cls.activation_data,
                                             {"v_init": -80, "t_stop": 40,
                                              "electrode_current": False,
                                              "chord_conductance":False,
                                              "normalization": "to_one"},
                                             1,
                                             "ActvationTestBKCa10uMCa",
                                             save_figures=True)
        cls.act_results = cls.test.run_model(cls.model,
                                             cls.test.stimulus_list,
                                             cls.test.v_init,
                                             cls.test.t_stop,
                                             cls.power,
                                             cls.test.t_mes,
                                             cls.test.chord_conductance,
                                             cls.test.electrode_current,
                                             "to_one")
       

    def test_summarize(self):
        self.score = self.test.judge(self.model)
        self.score.summarize()

    def test_run_model_keys(self):
        self.assertEqual(list(self.act_results.keys()),
                         self.test.stimulus_list)

    def test_run_model_values(self):
        values = np.array(list(self.act_results.values()))
        is_all_less_1 = np.all((values<=1))
        self.assertTrue(is_all_less_1)


class TestBK_100_uM_Ca(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.model = ModelCaConcClamp(channel_loc,
                                     ["bk"],
                                     ["k", "ca"],
                                     gbar_values={"bk": 1},
                                     external_conc={"k": 5},
                                     internal_conc={"ca": 100e-3,
                                                    "k": 120},
                                     temp=23, v_rest=-100,
                                     ljp=0)
        activation_data = np.loadtxt(activation_loc_BK_100_uM_Ca, skiprows=1,
                                     delimiter=",")
        cls.activation_data = {}
        for val in activation_data:
            cls.activation_data[val[0]] = val[1:].tolist()
        cls.power = 1
        cls.test = ActivationSteadyStateTest(cls.activation_data,
                                             {"v_init": -80, "t_stop": 40,
                                              "electrode_current": False,
                                              "chord_conductance":False,
                                              "normalization": "to_one"},
                                             1,
                                             "ActvationTestBKCa100uMCa",
                                             save_figures=True)
        cls.act_results = cls.test.run_model(cls.model,
                                             cls.test.stimulus_list,
                                             cls.test.v_init,
                                             cls.test.t_stop,
                                             cls.power,
                                             cls.test.t_mes,
                                             cls.test.chord_conductance,
                                             cls.test.electrode_current,
                                             "to_one")
       

    def test_summarize(self):
        self.score = self.test.judge(self.model)
        self.score.summarize()

    def test_run_model_keys(self):
        self.assertEqual(list(self.act_results.keys()),
                         self.test.stimulus_list)

    def test_run_model_values(self):
        values = np.array(list(self.act_results.values()))
        is_all_less_1 = np.all((values<=1))
        self.assertTrue(is_all_less_1)


        
if __name__ == "__main__":
    unittest.main()
