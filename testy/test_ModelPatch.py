import os
import pkg_resources
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

if __name__ == "__main__":
    unittest.main()
