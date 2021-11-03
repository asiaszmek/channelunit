import os
import unittest
from channelunit import ModelPatch

my_loc = os.path.dirname(os.path.abspath(__file__))
channel_loc = os.path.join(my_loc, "..", "demo_CA1", "ion_channels")


class TestModelPatch(unittest.TestCase):
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

class TestCapabilites(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.model = ModelPatch(channel_loc, "na3",
                                gbar_name="gbar", E_rev=40)
        cls.stim_levels_act = [-50, -40, -30, -20, -10, 0]
        cls.stim_levels_inact = [-105, -95, -85, -75, -65, -55, -45]
        cls.activation_cc = cls.model.get_activation_steady_state(
            cls.stim_levels_act, -90, 800, chord_conductance=True)
        cls.activation = cls.model.get_activation_steady_state(
            cls.stim_levels_act, -90, 800, chord_conductance=False)
        cls.inactivation_cc = cls.model.get_inactivation_steady_state(
            cls.stim_levels_inact, -5, 20, chord_conductance=True)
        cls.inactivation = cls.model.get_inactivation_steady_state(
            cls.stim_levels_inact, -5, 20, chord_conductance=False)
        
    def test_keys_activation(self):
        self.assertEqual(self.stim_levels_act,
                         list(self.activation.keys()))

    def test_keys_activation_cc(self):
        self.assertEqual(self.stim_levels_act,
                         list(self.activation_cc.keys()))

    def test_keys_inactivation(self):
        self.assertEqual(self.stim_levels_inact,
                         list(self.inactivation.keys()))

    def test_keys_inactivation_cc(self):
        self.assertEqual(self.stim_levels_inact,
                         list(self.inactivation_cc.keys()))


if __name__ == "__main__":
    unittest.main()
