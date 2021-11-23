import os
import unittest

import numpy as np

from channelunit.tests import SteadyStateTest
from channelunit.tests import InactivationSteadyStateTest
from channelunit.tests import ActivationSteadyStateTest
from channelunit import ModelPatch


my_loc = os.path.dirname(os.path.abspath(__file__))
channel_loc = os.path.join(my_loc, "..", "demo_CA1", "ion_channels")
activation_loc = os.path.join(my_loc, "..", "demo_CA1", "data",
                              "I_Na_activation.csv")
inactivation_loc = os.path.join(my_loc, "..", "demo_CA1", "data",
                                "I_Na_inactivation.csv")


class TestSteadyState(unittest.TestCase):
    def test_make_numeric_1(self):
        out = SteadyStateTest._get_numerical("5")
        self.assertEqual(out, 5)

    def test_make_numeric_2(self):
        out = SteadyStateTest._get_numerical(5)
        self.assertEqual(out, 5)

    def test_make_numeric_3(self):
        out = SteadyStateTest._get_numerical(5.)
        self.assertEqual(out, 5.)

    def test_make_numeric_4(self):
        self.assertRaises(ValueError, SteadyStateTest._get_numerical, "a")

    def test_make_numeric_5(self):
        self.assertRaises(ValueError, SteadyStateTest._get_numerical, [5, 7])

    def test_make_numeric_6(self):
        self.assertRaises(ValueError, SteadyStateTest._get_numerical, ["a"])
    
    def test_make_numeric_7(self):
        out = SteadyStateTest._get_numerical([5])
        self.assertEqual(out, 5)

    def test_make_numeric_8(self):
        out = SteadyStateTest._get_numerical([5.])
        self.assertEqual(out, 5.)

    def test_extract_simulation(self):
        obs = {"-30": 3, -50: 6, -10.: 1}
        out = SteadyStateTest.extract_stimulation(obs)
        expected = [-50, -30, -10.]
        self.assertEqual(expected, out)

    def test_extract_simulation_1(self):
        obs = {"-30": 3, -50:6, -10.:1, "a":2}
        self.assertRaises(ValueError, SteadyStateTest.extract_stimulation,
                          obs)

    def test_format_data_1(self):
        obs = {"-10": "3", "-20": ["2", "0.5"]}
        out = SteadyStateTest.format_data(obs)
        self.assertEqual(out, {-10:[3, 0.03], -20:[2, 0.5]})


class TestActivationSteadyState(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        E_rev = 8.314*(273.15+22)/96485*np.log(110/15)
        cls.model = ModelPatch(channel_loc, "na3", "na", external_conc=110, temp=22,
                               liquid_junction_pot=0)
        activation_data = np.loadtxt(activation_loc, skiprows=1,
                                     delimiter=",")
        cls.power = 2
        cls.activation_data = dict(val.tolist() for val in activation_data)
        cls.test = ActivationSteadyStateTest(cls.activation_data,
                                             {"v_init": -90, "t_stop": 400,
                                              "chord_conductance":True}, 2,
                                             "ActvationSSTest",
                                             save_figures=True)



    def test_summarize(self):
        self.score = self.test.judge(self.model)
        self.score.summarize()
        
    def test_run_model(self):
        out = self.test.run_model(self.model, self.test.stimulus_list,
                                  self.test.v_init, self.test.t_stop,
                                  self.power,
                                  self.test.chord_conductance)
        self.assertEqual(list(out.keys()), self.test.stimulus_list)


class TestInactivationSteadyState(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        E_rev = 8.314*(273.15+22)/96485*np.log(110/15)
        cls.model = ModelPatch(channel_loc, "na3", "na", external_conc=110,
                               temp=22,
                               liquid_junction_pot=0)
        inactivation_data = np.loadtxt(inactivation_loc, skiprows=1,
                                     delimiter=",")
        cls.inactivation_data = dict(val.tolist() for val in inactivation_data)
        cls.power = 1
        cls.test = InactivationSteadyStateTest(cls.inactivation_data,
                                               {"v_test": -5, "t_test": 100,
                                                "chord_conductance":True}, 1,
                                               "InactvationSSTest",
                                               save_figures=True)



    def test_summarize(self):
        self.score = self.test.judge(self.model)
        self.score.summarize()
        
    def test_run_model(self):
        out = self.test.run_model(self.model, self.test.stimulus_list,
                                  self.test.v_test, self.test.t_test, self.power,
                                  self.test.chord_conductance)
        self.assertEqual(list(out.keys()), self.test.stimulus_list)

if __name__ == "__main__":
    unittest.main()
