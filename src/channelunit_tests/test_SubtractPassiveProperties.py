import unittest
import numpy as np

from channelunit.base_classes import MembranePatch


class TestSubtractPassiveProperties(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.modelljp2 = MembranePatch(temp=22, Rm=20000, cm=2,
                                      v_rest=-65, ljp=10,
                                      cvode=True, sim_dt=0.0001)
        cls.modelljp2.patch.L = 10
        cls.modelljp2.patch.diam= 10
        cls.dur1 = 20
        cls.dur2 = 20
        cls.delay = cls.dur2
        cls.modelljp2.set_vclamp(cls.dur1, 10, cls.dur2, 90, True)
        cls.pulse = 20
        t_stop = cls.dur1 + 6*cls.dur2 + 5*cls.delay
        cls.current = cls.modelljp2.run(t_stop)


    def test_current_subtraction(self):
        dt = self.modelljp2.dt
        passive = (self.current[int(self.dur1/dt):
                               int((self.dur1+self.dur2)/dt)].copy()-
                   self.current[:int(self.dur1/dt)].copy())
        expected = self.modelljp2.curr_stim_response(self.current,
                                                     self.dur1,
                                                     self.dur2, dt, True)
        self.assertTrue(np.allclose(passive[4:], expected[4:]))

    def test_pulse_sum(self):
        dt = self.modelljp2.dt
        length = int(self.dur2/self.modelljp2.dt)
        p_sum = np.zeros((length,))
        t_start = self.dur1 + 2*self.dur2 
        for i in range(4):
            basal = self.current[int(t_start/dt):
                                  int((t_start+self.dur2)/dt)]
            p_sum += (self.current[int((t_start+self.dur2)/dt):
                                  int((t_start+2*self.dur2)/dt)]-basal)
            t_start += 2*self.dur2
        expected = self.modelljp2.curr_leak_amp(self.current, self.dur1,
                                                self.dur2, dt)
        self.assertTrue(np.allclose(p_sum, expected))

    def test_leak_substraction(self):
        dt = self.modelljp2.dt
        expected = (self.current[int(self.dur1/dt):
                                int((self.dur1+self.dur2)/dt)].copy()-
                    self.current[:int(self.dur1/dt)].copy())
        automatic_leak_subtraction = self.modelljp2.curr_leak_amp(self.current,
                                                                  self.dur1,
                                                                  self.dur2,
                                                                  dt)
        
        difference = abs(expected - automatic_leak_subtraction)
        self.assertTrue(np.all(difference[4:] < 0.03))

    def test_sum_amp(self):
        out_1 = (self.modelljp2.vclamp.amp5 + self.modelljp2.vclamp.amp7
                 + self.modelljp2.vclamp.amp9 + self.modelljp2.vclamp.amp11)
        out_2 = 4*self.pulse
        self.assertEqual(out_1-4*self.modelljp2.vclamp.amp10, out_2)

    def test_sum_amp_2(self):
        out_1 = (self.modelljp2.vclamp.amp5 + self.modelljp2.vclamp.amp7
                 + self.modelljp2.vclamp.amp9 + self.modelljp2.vclamp.amp11)
        out_2 = self.modelljp2.vclamp.amp2 - self.modelljp2.vclamp.amp1
        self.assertEqual(out_1-4*self.modelljp2.vclamp.amp10, out_2)


if __name__ == "__main__":
    unittest.main()

