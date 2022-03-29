import numpy
from sciunit import Score
from sciunit.errors import InvalidScoreError

class ZScore_SteadyStateCurves(Score):
    """
    ZScore of all the points of the activation/inactivation curve.
    """
    def __init__(self, score, related_data={}):
        if not isinstance(score, Exception) and not isinstance(score, float):
            raise InvalidScoreError("Score must be a float.")
        else:
            super(ZScore_SteadyStateCurves, self).__init__(score,
                                                          related_data=related_data)

    @classmethod
    def compute(cls, observation, prediction):
        errors = {}
        for key in prediction.keys():
            error = abs(prediction[key] - observation[key][0])/observation[key][1]
            errors[key] = error
        print(errors)
        print(numpy.nanmean(list(errors.values())))
        return numpy.nanmean(list(errors.values())), errors

    def __str__(self):
        return "Avg ZScore=%f" % self.score


class ZScore_BothSteadyStateCurves(Score):
    """
    ZScore of all the points of the activation/inactivation curve.
    """
    def __init__(self, score, related_data={}):
        if not isinstance(score, Exception) and not isinstance(score, float):
            raise InvalidScoreError("Score must be a float.")
        else:
            super(ZScore_BothSteadyStateCurves, self).__init__(score,
                                                               related_data=related_data)

    @classmethod
    def compute(cls, observation, prediction):
        errors = {}
        assert sorted(prediction.keys()) == sorted(observation.keys())
        assert sorted(prediction.keys()) == ["Activation", "Inactivation"]
        for exp in ["Activation", "Inactivation"]:
            for key in prediction[exp].keys():
                error = abs(prediction[exp][key] - observation[exp][key][0])/observation[exp][key][1]
                new_key = "%s_%s" % (exp, key)
                errors[new_key] = error
        print(errors)
        print(numpy.nanmean(list(errors.values())))
        return numpy.nanmean(list(errors.values())), errors

    def __str__(self):
        return "Avg ZScore=%f" % self.score
