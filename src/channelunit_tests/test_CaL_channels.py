import os
import unittest

import numpy as np

from channelunit.tests import ActivationSteadyStateTest
from channelunit import ModelWholeCellPatchCaSingleChan
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

class TestCaLChannelsLowBariumba(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.modelBa20_H = ModelWholeCellPatchCaSingleChan(channel_loc,
                                                               "calHGHK",
                                                               "ba",
                                                               external_conc=20,
                                                               temp=22,
                                                               ljp=0)
        activation_data = np.loadtxt(activation_loc_Cal_20_Ba, skiprows=1,
                                     delimiter=",")
        cls.power = 1
        cls.activation_data = dict(val.tolist() for val in activation_data)
        cls.test_Ba20 = ActivationSteadyStateTest(cls.activation_data,
                                                  {"v_init": -90, "t_stop": 70,
                                                   "electrode_current": False,
                                                   "chord_conductance":False,
                                                   "normalization": "to_one"}, 1,
                                                  "ActvationSSTestLowba",
                                                  save_figures=True)
        cls.act_results = cls.test_Ba20.run_model(cls.modelBa20_H,
                                                  cls.test_Ba20.stimulus_list,
                                                  cls.test_Ba20.v_init,
                                                  cls.test_Ba20.t_stop,
                                                  cls.power,
                                                  cls.test_Ba20.chord_conductance,
                                                  cls.test_Ba20.electrode_current,
                                                  "to_one")
       

    def test_summarize_H(self):
        self.score = self.test_Ba20.judge(self.modelBa20_H)
        self.score.summarize()

    def test_run_model_H_keys(self):
        self.assertEqual(list(self.act_results.keys()),
                         self.test_Ba20.stimulus_list)

    def test_run_model_H_values(self):
        values = np.array(list(self.act_results.values()))
        is_all_less_1 = np.all((values<=1))
        self.assertTrue(is_all_less_1)

        
class TestCaLChannelsLowBariumBa(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        
        cls.modelBa20_HH = ModelWholeCellPatchCaSingleChan(channel_loc,
                                                      "CalHGHK",
                                                      "Ba",
                                                      external_conc=20,
                                                      temp=22,
                                                      ljp=0)

        activation_data = np.loadtxt(activation_loc_Cal_20_Ba, skiprows=1,
                                     delimiter=",")
        cls.power = 1
        cls.activation_data = dict(val.tolist() for val in activation_data)
        cls.test_Ba20 = ActivationSteadyStateTest(cls.activation_data,
                                                  {"v_init": -90, "t_stop": 70,
                                                   "electrode_current": False,
                                                   "chord_conductance": False,
                                                   "normalization": "to_one"}, 1,
                                                  "ActvationSSTestLowBa",
                                                  save_figures=True)
        cls.act_results = cls.test_Ba20.run_model(cls.modelBa20_HH,
                                                  cls.test_Ba20.stimulus_list,
                                                  cls.test_Ba20.v_init,
                                                  cls.test_Ba20.t_stop,
                                                  cls.power,
                                                  cls.test_Ba20.chord_conductance,
                                                  cls.test_Ba20.electrode_current,
                                                  "save_sign")


    def test_summarize_H(self):
        self.score = self.test_Ba20.judge(self.modelBa20_HH)
        self.score.summarize()

    def test_run_model_H_keys(self):
        self.assertEqual(list(self.act_results.keys()),
                         self.test_Ba20.stimulus_list)

    def test_run_model_H_values(self):
        values = np.array(list(self.act_results.values()))
        is_all_less_1 = np.all((values<=0))
        self.assertTrue(is_all_less_1)



class TestCaLChannelsLowCalciumca(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.modelCa20_H = ModelWholeCellPatchCaSingleChan(channel_loc,
                                                     "calHGHK",
                                                     "ca",
                                                     external_conc=1.5,
                                                     temp=22,
                                                     ljp=0)
        activation_data = np.loadtxt(activation_loc_Cal12, skiprows=1,
                                     delimiter=",")
        cls.power = 1
        cls.activation_data = dict()
        for val in activation_data:
            cls.activation_data[val[0]] = val[1:].tolist()
        cls.test_Ca20 = ActivationSteadyStateTest(cls.activation_data,
                                                  {"v_init": -90, "t_stop": 70,
                                                   "electrode_current": True,
                                                   "chord_conductance":False,
                                                   "normalization":
                                                   "save_sign"}, 1,
                                                  "ActvationSSTestLowca",
                                                  save_figures=True)
        cls.act_results = cls.test_Ca20.run_model(cls.modelCa20_H,
                                                  cls.test_Ca20.stimulus_list,
                                                  cls.test_Ca20.v_init,
                                                  cls.test_Ca20.t_stop,
                                                  cls.power,
                                                  cls.test_Ca20.chord_conductance,
                                                  cls.test_Ca20.electrode_current,
                                                  "save_sign")
       

    def test_summarize_H(self):
        self.score = self.test_Ca20.judge(self.modelCa20_H)
        self.score.summarize()

    def test_run_model_H_keys(self):
        self.assertEqual(list(self.act_results.keys()),
                         self.test_Ca20.stimulus_list)

    def test_run_model_H_values(self):
        values = np.array(list(self.act_results.values()))
        is_all_less_1 = np.all((values<=1))
        self.assertTrue(is_all_less_1)

        
class TestCaLChannelsLowCalciumCa(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        
        cls.modelCa20_HH = ModelWholeCellPatchCaSingleChan(channel_loc,
                                                                "CalHGHK",
                                                                "Ca",
                                                                external_conc=1.5,
                                                                temp=22,
                                                                ljp=0)

        activation_data = np.loadtxt(activation_loc_Cal12, skiprows=1,
                                     delimiter=",")
        cls.power = 1
        cls.activation_data = dict()
        for val in activation_data:
            cls.activation_data[val[0]] = val[1:].tolist()
        cls.test_Ca20 = ActivationSteadyStateTest(cls.activation_data,
                                                {"v_init": -90, "t_stop": 70,
                                                 "electrode_current": False,
                                                 "chord_conductance":False,
                                                 "normalization": "save_sign"}, 1,
                                                "ActvationSSTestLowCa",
                                                save_figures=True)
        cls.act_results = cls.test_Ca20.run_model(cls.modelCa20_HH,
                                                  cls.test_Ca20.stimulus_list,
                                                  cls.test_Ca20.v_init,
                                                  cls.test_Ca20.t_stop,
                                                  cls.power,
                                                  cls.test_Ca20.chord_conductance,
                                                  cls.test_Ca20.electrode_current,
                                                  "save_sign")


    def test_summarize_H(self):
        self.score = self.test_Ca20.judge(self.modelCa20_HH)
        self.score.summarize()

    def test_run_model_H_keys(self):
        self.assertEqual(list(self.act_results.keys()),
                         self.test_Ca20.stimulus_list)

    def test_run_model_H_values(self):
        values = np.array(list(self.act_results.values()))
        
        is_all_less_1 = np.all((values<=1))

        self.assertTrue(is_all_less_1)


if __name__ == "__main__":
    unittest.main()
