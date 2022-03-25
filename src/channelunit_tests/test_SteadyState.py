import os
import unittest

import numpy as np

from channelunit.tests import BaseSteadyStateTest
from channelunit.tests import InactivationSteadyStateTest
from channelunit.tests import ActivationSteadyStateTest
from channelunit import ModelWholeCellPatchSingleChan
from channelunit import data_path


channel_loc = os.path.join(data_path, "ion_channels")
activation_loc = os.path.join(data_path, "data",
                              "I_Na_activation.csv")
inactivation_loc = os.path.join(data_path, "data",
                                "I_Na_inactivation.csv")


class TestBaseSteadyState(unittest.TestCase):
    def test_make_numeric_1(self):
        out = BaseSteadyStateTest._get_numerical("5")
        self.assertEqual(out, 5)

    def test_make_numeric_2(self):
        out = BaseSteadyStateTest._get_numerical(5)
        self.assertEqual(out, 5)

    def test_make_numeric_3(self):
        out = SteadyStateTest._get_numerical(5.)
        self.assertEqual(out, 5.)

    def test_make_numeric_4(self):
        self.assertRaises(ValueError, BaseSteadyStateTest._get_numerical, "a")

    def test_make_numeric_5(self):
        self.assertRaises(ValueError, BaseSteadyStateTest._get_numerical, [5, 7])

    def test_make_numeric_6(self):
        self.assertRaises(ValueError, BaseSteadyStateTest._get_numerical, ["a"])
    
    def test_make_numeric_7(self):
        out = BaseSteadyStateTest._get_numerical([5])
        self.assertEqual(out, 5)

    def test_make_numeric_8(self):
        out = BaseSteadyStateTest._get_numerical([5.])
        self.assertEqual(out, 5.)

    def test_extract_simulation(self):
        obs = {"-30": 3, -50: 6, -10.: 1}
        out = BaseSteadyStateTest.extract_stimulation(obs)
        expected = [-50, -30, -10.]
        self.assertEqual(expected, out)

    def test_extract_simulation_1(self):
        obs = {"-30": 3, -50:6, -10.:1, "a":2}
        self.assertRaises(ValueError, BaseSteadyStateTest.extract_stimulation,
                          obs)

    def test_format_data_1(self):
        obs = {"-10": "3", "-20": ["2", "0.5"]}
        out = BaseSteadyStateTest.format_data(obs)
        self.assertEqual(out, {-10:[3, 0.03], -20:[2, 0.5]})


class TestActivationSteadyState(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        E_rev = 8.314*(273.15+22)/96485*np.log(110/15)
        cls.model = ModelWholeCellPatchSingleChan(channel_loc, "na3", "na",
                                     external_conc=110, temp=22,
                                     ljp=0)
        activation_data = np.loadtxt(activation_loc, skiprows=1,
                                     delimiter=",")
        cls.power = 1
        cls.activation_data = dict(val.tolist() for val in activation_data)
        cls.test_no_ls = ActivationSteadyStateTest(cls.activation_data,
                                             {"v_init": -90, "t_stop": 200,
                                              "electrode_current":False,
                                              "chord_conductance":True,
                                              "normalization": "to_one"}, 1,
                                                   "ActvationSSTest",
                                                   save_figures=True)
        cls.test_ls = ActivationSteadyStateTest(cls.activation_data,
                                                {"v_init": -90, "t_stop": 200,
                                                 "electrode_current":True,
                                                 "chord_conductance":True,
                                                 "normalization": "to_one"}, 1,
                                                "ActvationSSTest",
                                                save_figures=True)



    def test_summarize(self):
        self.score = self.test_ls.judge(self.model)
        self.score.summarize()


    def test_summarize_no_ls(self):
        self.score = self.test_no_ls.judge(self.model)
        self.score.summarize()

    def test_run_model_1_keys(self):
        out = self.test_no_ls.run_model(self.model,
                                        self.test_no_ls.stimulus_list,
                                        self.test_no_ls.v_init,
                                        self.test_no_ls.t_stop,
                                        self.power,
                                        self.test_no_ls.chord_conductance,
                                        self.test_no_ls.electrode_current,
                                        self.test_no_ls.normalization)
        self.assertEqual(list(out.keys()), self.test_no_ls.stimulus_list)

    def test_run_model_1_values(self):
        out = self.test_no_ls.run_model(self.model, self.test_no_ls.stimulus_list,
                                        self.test_no_ls.v_init,
                                        self.test_no_ls.t_stop,
                                        self.power,
                                        self.test_no_ls.chord_conductance,
                                        self.test_no_ls.electrode_current,
                                        self.test_no_ls.normalization)
        
        values = np.array(list(out.values()))
        is_all_less_1 = np.all((values<=1))
        self.assertTrue(is_all_less_1)

    def test_run_model_2_keys(self):
        out = self.test_ls.run_model(self.model, self.test_ls.stimulus_list,
                                     self.test_ls.v_init, self.test_ls.t_stop,
                                     self.power,
                                     self.test_ls.chord_conductance,
                                     self.test_ls.electrode_current,
                                     self.test_ls.normalization)
        self.assertEqual(list(out.keys()), self.test_ls.stimulus_list)

    def test_run_model_2_values(self):
        out = self.test_ls.run_model(self.model, self.test_ls.stimulus_list,
                                     self.test_ls.v_init,
                                     self.test_ls.t_stop,
                                     self.power,
                                     self.test_ls.chord_conductance,
                                     self.test_ls.electrode_current,
                                     self.test_ls.normalization)
        values = np.array(list(out.values()))
        is_all_less_1 = np.all((values<=1))
        self.assertTrue(is_all_less_1)


class TestInactivationSteadyState(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        E_rev = 8.314*(273.15+22)/96485*np.log(110/15)
        cls.model = ModelWholeCellPatchSingleChan(channel_loc, "na3", "na",
                                                  external_conc=110,
                                                  temp=22,
                                                  ljp=0)
        inactivation_data = np.loadtxt(inactivation_loc, skiprows=1,
                                       delimiter=",")
        cls.inactivation_data = dict(val.tolist() for val in inactivation_data)
        cls.power = 1
        cls.test = InactivationSteadyStateTest(cls.inactivation_data,
                                               {"v_test": -5, "t_test": 10,
                                                "chord_conductance":True,
                                                "electrode_current":True,
                                                "normalization": "to_one"}, 1,
                                               "InactivationSSTest",
                                               save_figures=True)



    def test_summarize(self):
        self.score = self.test.judge(self.model)
        self.score.summarize()
        
    def test_run_model(self):
        out = self.test.run_model(self.model, self.test.stimulus_list,
                                  self.test.v_test, self.test.t_test,
                                  self.power,
                                  self.test.chord_conductance,
                                  self.test.electrode_current,
                                  self.test.normalization)
        self.assertEqual(list(out.keys()), self.test.stimulus_list)

if __name__ == "__main__":
    unittest.main()
