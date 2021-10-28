import unittest
from channelunit.scores import ZScore_SteadyStateCurves as ZScore_SS

class Test_ZScore_SteadyStateCurves(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        observation = {-70: [1, 0.02],
                       -60: [2, 1]}
        prediction = {-70: 0,
                      -60: 5}
        cls.score, cls.errors = ZScore_SS.compute(observation, prediction)

    def test_calc(self):
        out = (1/0.02+3)/2
        self.assertEqual(self.score, out)

    def test_errors(self):
        errors = {-70: 1/0.02, -60: 3}
        self.assertEqual(self.errors, errors)
        
if __name__ == "__main__":
    unittest.main()
