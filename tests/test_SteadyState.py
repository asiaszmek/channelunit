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


if __name__ == "__main__":
    unittest.main()
