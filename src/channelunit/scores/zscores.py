import numpy
from sciunit import Score

class ZScore_SteadyStateCurves(Score):
    """
    ZScore of all the points of the activation/inactivation curve.
    """
    def __init__(self, score, related_data={}):
        if not isinstance(score, Exception) and not isinstance(score, float):
            raise InvalidScoreError("Score must be a float.")
        else:
            super(ZScore_SteadyStateCurves,self).__init__(score,
                                                          related_data=related_data)

    @classmethod
    def compute(cls, observation, prediction):
        errors = {}
        for key in prediction.keys():
            error = abs(prediction[key] - observation[key][0])/observation[key][1]
            errors[key] = error
        return numpy.nanmean(list(errors.values())), errors

    def __str__(self):
        return "Avg ZScore=%f" % self.score
