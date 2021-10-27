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
if __name__ == "__main__":
    unittest.main()
