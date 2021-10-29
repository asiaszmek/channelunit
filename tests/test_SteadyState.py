import unittest
from channelunit.tests import SteadyStateTest

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

    
if __name__ == "__main__":
    unittest.main()
