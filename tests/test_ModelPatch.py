import os
import unittest

import numpy as np

from channelunit import ModelPatch

my_loc = os.path.dirname(os.path.abspath(__file__))
channel_loc = os.path.join(my_loc, "..", "demo_CA1", "ion_channels")


class TestModelPatch(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.modelJ = ModelPatch(channel_loc, "nap", gbar_name="gnabar")
        cls.modelNJ = ModelPatch(channel_loc, "nap", gbar_name="gnabar",
                                 liquid_junction_pot=0)
        cls.modelJ.set_vclamp(10, 10, 100, 100)
        cls.modelNJ.set_vclamp(10, 10, 100, 100)

    def test_reading_in_easy(self):
        out = ModelPatch(channel_loc, "na3")
        self.assertEqual(50, out.E_rev)

    def test_reading_in_no_gbar(self):
        self.assertRaises(SystemExit,  ModelPatch,
                          channel_loc, "na3", gbar_name="gbar1")

    def test_reading_in_provide_E_rev(self):
        out = ModelPatch(channel_loc, "na3", E_rev=40)
        self.assertEqual(40, out.E_rev)

    def test_no_E_rev_name(self):
        self.assertRaises(SystemExit,  ModelPatch,
                          channel_loc, "hd")

    def test_no_E_rev_name_provide_E_rev(self):
        out = ModelPatch(channel_loc, "hd", E_rev=-30)
        self.assertEqual(-30, out.E_rev)

    def test_setup_gbar(self):
        out = ModelPatch(channel_loc, "nax", E_rev=-30)
        self.assertEqual(0.001,
                         out.soma.psection()["density_mechs"]["nax"]["gbar"][0])

    def test_setup_gbar_custom(self):
        out = ModelPatch(channel_loc, "nap", gbar_name="gnabar")
        self.assertEqual(0.001,
                         out.soma.psection()["density_mechs"]["nap"]["gnabar"][0])

    def test_get_E_rev_value(self):
        out = ModelPatch(channel_loc, "nap", gbar_name="gnabar")
        val = out.get_E_rev_value()
        self.assertEqual(50, val)
        
    
    def test_get_E_rev_name(self):
        out = ModelPatch(channel_loc, "nap", gbar_name="gnabar")
        name = out.get_E_rev_name()
        self.assertEqual("ena", name)

    def test_set_vclamp_junction_amp1(self):
        self.assertEqual(self.modelJ.vclamp.amp1,
                         10-self.modelJ.junction)

    def test_set_vclamp_junction_amp2(self):
        self.assertEqual(self.modelJ.vclamp.amp2,
                         100-self.modelJ.junction)

    def test_set_vclamp_nojunction_amp1(self):
        self.assertEqual(self.modelNJ.vclamp.amp1,
                         10)

    def test_set_vclamp_nojunction_amp2(self):
        self.assertEqual(self.modelNJ.vclamp.amp2,
                         100)

    def test_set_vclamp_junction_dur1(self):
        self.assertEqual(self.modelJ.vclamp.dur1,
                         10)

    def test_set_vclamp_junction_dur2(self):
        self.assertEqual(self.modelJ.vclamp.dur2,
                         100)

    def test_set_vclamp_nojunction_dur1(self):
        self.assertEqual(self.modelNJ.vclamp.dur1,
                         10)

    def test_set_vclamp_nojunction_dur2(self):
        self.assertEqual(self.modelNJ.vclamp.dur2,
                         100)

    def test_extract_current(self):
        I = 2
        out = self.modelJ.extract_current(I, True)
        expected = I/(self.modelJ.vclamp.amp2 - self.modelJ.E_rev)
        self.assertEqual(out, expected)

    def test_extract_current(self):
        I = -22
        out = self.modelJ.extract_current(I, False)
        self.assertEqual(out, 22)

    def test_normalize_to_one(self):
        dic = {1:1, 2:2}
        out = self.modelJ.normalize_to_one(dic)
        self.assertEqual({1:0.5, 2:1.}, out)
        
class TestCapabilites(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.modelY = ModelPatch(channel_loc, "na3",
                                gbar_name="gbar", E_rev=40, cvode=True)
        cls.stim_levels_act = [-50, -40, -30, -20, -10, 0]
        cls.stim_levels_inact = [-105, -95, -85, -75, -65, -55, -45]
        cls.activationY_cc = cls.modelY.get_activation_steady_state(
            cls.stim_levels_act, -90, 800, chord_conductance=True)
        cls.activationY = cls.modelY.get_activation_steady_state(
            cls.stim_levels_act, -90, 800, chord_conductance=False)
        cls.inactivationY_cc = cls.modelY.get_inactivation_steady_state(
            cls.stim_levels_inact, -5, 20, chord_conductance=True)
        cls.inactivationY = cls.modelY.get_inactivation_steady_state(
            cls.stim_levels_inact, -5, 20, chord_conductance=False)

        cls.modelN = ModelPatch(channel_loc, "na3",
                                gbar_name="gbar", E_rev=40, cvode=False)
        cls.stim_levels_act = [-50, -40, -30, -20, -10, 0]
        cls.stim_levels_inact = [-105, -95, -85, -75, -65, -55, -45]
        cls.activationN_cc = cls.modelN.get_activation_steady_state(
            cls.stim_levels_act, -90, 800, chord_conductance=True)
        cls.activationN = cls.modelN.get_activation_steady_state(
            cls.stim_levels_act, -90, 800, chord_conductance=False)
        cls.inactivationN_cc = cls.modelN.get_inactivation_steady_state(
            cls.stim_levels_inact, -5, 20, chord_conductance=True)
        cls.inactivationN = cls.modelN.get_inactivation_steady_state(
            cls.stim_levels_inact, -5, 20, chord_conductance=False)

    def test_keys_activationY(self):
        self.assertEqual(self.stim_levels_act,
                         list(self.activationY.keys()))

    def test_keys_activationY_cc(self):
        self.assertEqual(self.stim_levels_act,
                         list(self.activationY_cc.keys()))

    def test_keys_inactivationY(self):
        self.assertEqual(self.stim_levels_inact,
                         list(self.inactivationY.keys()))

    def test_keys_inactivationY_cc(self):
        self.assertEqual(self.stim_levels_inact,
                         list(self.inactivationY_cc.keys()))

    def test_keys_activationN(self):
        self.assertEqual(self.stim_levels_act,
                         list(self.activationN.keys()))

    def test_keys_activationN_cc(self):
        self.assertEqual(self.stim_levels_act,
                         list(self.activationN_cc.keys()))

    def test_keys_inactivationN(self):
        self.assertEqual(self.stim_levels_inact,
                         list(self.inactivationN.keys()))

    def test_keys_inactivationN_cc(self):
        self.assertEqual(self.stim_levels_inact,
                         list(self.inactivationN_cc.keys()))
        

if __name__ == "__main__":
    unittest.main()
