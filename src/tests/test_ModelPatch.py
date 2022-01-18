import os
import unittest

import numpy as np

from channelunit.model_patch import ModelPatch


class TestSubtractPassiveProperties(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.modelljp1 = ModelPatch(liquid_junction_pot=10)
        cls.modelljp2 = ModelPatch(liquid_junction_pot=10)
        cls.model =  ModelPatch()
        cls.modelljp1.set_vclamp(10, 10, 100, 100, False)
        cls.dur1 = 10
        cls.dur2 = 20
        cls.delay = 30
        cls.modelljp2.set_vclamp(cls.dur1, 10, cls.dur2, 90, True, cls.delay)
        cls.pulse = 20

    def test_set_vclamp_junction_amp1(self):
        self.assertEqual(self.modelljp1.vclamp.amp1,
                         10-self.modelljp1.junction)

    def test_set_vclamp_junction_amp2(self):
        self.assertEqual(self.modelljp1.vclamp.amp2,
                         100-self.modelljp1.junction)

          
    def test_set_vclamp_junction_dur1(self):
        self.assertEqual(self.modelljp1.vclamp.dur1, 10)

    def test_set_vclamp_junction_dur2(self):
        self.assertEqual(self.modelljp1.vclamp.dur2,
                         100)

    def test_set_vclamp_ls_junction_amp1(self):
        self.assertEqual(self.modelljp2.vclamp.amp1,
                         10-self.modelljp2.junction)

    def test_set_vclamp_ls_junction_amp2(self):
        self.assertEqual(self.modelljp2.vclamp.amp2,
                         90-self.modelljp2.junction)
    
    def test_set_vclamp_ls_junction_dur1(self):
        self.assertEqual(self.modelljp2.vclamp.dur1,
                         10)

    def test_set_vclamp_ls_junction_dur2(self):
        self.assertEqual(self.modelljp2.vclamp.dur2,
                         20)

    def test_set_vclamp_ls_junction_amp3(self):
        self.assertEqual(self.modelljp2.vclamp.amp3,
                         10-self.modelljp2.junction)

    def test_set_vclamp_ls_junction_dur3(self):
        self.assertEqual(self.modelljp2.vclamp.dur3,
                         30)

    def test_set_vclamp_ls_junction_amp4(self):
        self.assertEqual(self.modelljp2.vclamp.amp4,
                         10-2*self.modelljp2.junction-
                         self.modelljp2.v_low-self.pulse)

    def test_set_vclamp_ls_junction_dur4(self):
        self.assertEqual(self.modelljp2.vclamp.dur4,
                         30)

    def test_set_vclamp_ls_junction_amp5(self):
        self.assertEqual(self.modelljp2.vclamp.amp5,
                         10-2*self.modelljp2.junction-self.modelljp2.v_low)

    def test_set_vclamp_ls_junction_dur5(self):
        self.assertEqual(self.modelljp2.vclamp.dur5,
                         20)

    def test_set_vclamp_ls_junction_amp6(self):

        self.assertEqual(self.modelljp2.vclamp.amp6,
                         10-2*self.modelljp2.junction-self.modelljp2.v_low-
                         self.pulse)

    def test_set_vclamp_ls_junction_dur6(self):
        self.assertEqual(self.modelljp2.vclamp.dur6,
                         30)

    def test_set_vclamp_ls_junction_amp7(self):
        self.assertEqual(self.modelljp2.vclamp.amp7,
                         10-2*self.modelljp2.junction-self.modelljp2.v_low)

    def test_set_vclamp_ls_junction_dur7(self):
        self.assertEqual(self.modelljp2.vclamp.dur7,
                         20)

    def test_set_vclamp_ls_junction_amp8(self):
        self.assertEqual(self.modelljp2.vclamp.amp8,
                         10-2*self.modelljp2.junction-self.modelljp2.v_low
                         -self.pulse)

    def test_set_vclamp_ls_junction_dur8(self):
        self.assertEqual(self.modelljp2.vclamp.dur8,
                         30)

    def test_set_vclamp_ls_junction_amp9(self):
        self.assertEqual(self.modelljp2.vclamp.amp9,
                         10-2*self.modelljp2.junction-self.modelljp2.v_low)

    def test_set_vclamp_ls_junction_dur9(self):
        self.assertEqual(self.modelljp2.vclamp.dur9,
                         20)
    def test_set_vclamp_ls_junction_amp10(self):
        self.assertEqual(self.modelljp2.vclamp.amp10,
                         10-2*self.modelljp2.junction-self.modelljp2.v_low
                         -self.pulse)

    def test_set_vclamp_ls_junction_dur10(self):
        self.assertEqual(self.modelljp2.vclamp.dur10,
                         30)

    def test_set_vclamp_ls_junction_amp11(self):
        self.assertEqual(self.modelljp2.vclamp.amp11,
                         10-2*self.modelljp2.junction-self.modelljp2.v_low)

    def test_set_vclamp_ls_junction_dur11(self):
        self.assertEqual(self.modelljp2.vclamp.dur11,
                         20)

    def test_set_vclamp_ls_junction_amp12(self):
        self.assertEqual(self.modelljp2.vclamp.amp12,
                         10-2*self.modelljp2.junction-self.modelljp2.v_low
                         -self.pulse)

    def test_set_vclamp_ls_junction_dur12(self):
        self.assertEqual(self.modelljp2.vclamp.dur12,
                         30)

    def test_pulse_height_1(self):
        self.assertEqual(self.modelljp2.vclamp.amp5,
                         self.modelljp2.vclamp.amp4+self.pulse)

    def test_pulse_height_2(self):
        self.assertEqual(self.modelljp2.vclamp.amp7,
                         self.modelljp2.vclamp.amp6+self.pulse)

    def test_pulse_height_3(self):
        self.assertEqual(self.modelljp2.vclamp.amp9, 
                         self.modelljp2.vclamp.amp8+self.pulse)
        
    def test_pulse_height_4(self):
        self.assertEqual(self.modelljp2.vclamp.amp11,
                         self.modelljp2.vclamp.amp10+self.pulse)

    def test_leak_subtraction(self):
        t_stop = self.dur1 + 6*self.dur2 +5*self.delay
        current = self.modelljp2.run(t_stop)
        dt = self.modelljp2.dt
        passive = current[int(self.dur1/dt)+1: int((self.dur1+self.dur2)/dt)]
        p_sum = np.zeros(passive.shape)
        t_start = self.dur1 + self.dur2 + 2*self.delay
        for i in range(4):
            p_sum += current[int(t_start/dt)+1: int((t_start+self.dur2)/dt)]
            t_start += self.dur2 + self.delay
        new_current = passive - p_sum
        expected = np.zeros(new_current.shape)
        print(new_current)
        self.assertTrue(np.allclose(new_current, expected))

if __name__ == "__main__":
    unittest.main()

