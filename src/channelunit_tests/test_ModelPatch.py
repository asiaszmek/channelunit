import os
import unittest

import numpy as np

from channelunit.model_patch import ModelPatch

class TestVclamp(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.modelljp1 = ModelPatch(ljp=10)
        cls.modelljp2 = ModelPatch(ljp=10, cm=2)
        cls.model =  ModelPatch()
        cls.modelljp1.set_vclamp(10, 10, 100, 100, False)
        cls.dur1 = 10
        cls.dur2 = 20
        cls.delay = 30
        cls.modelljp2.set_vclamp(cls.dur1, 10, cls.dur2, 90, True,
                                 cls.delay)
        cls.pulse = 20

    def test_init_default_1(self):
        self.assertEqual(self.modelljp1.temperature, 22)

    def test_init_default_pas_1(self):
        self.assertTrue(np.isclose(self.modelljp1.patch.g_pas,
                                   1/20000))

    def test_init_default_pas_2(self):
        self.assertTrue(np.isclose(self.modelljp1.patch.e_pas,
                                   -65))

    def test_init_default_cvode(self):
        self.assertTrue(self.modelljp1.cvode)

        
    def test_set_cm(self):
        cm = self.modelljp2.patch.cm
        self.assertEqual(cm, 2)

    def test_setter_cm(self):
        self.modelljp1.cm = 2
        self.assertEqual(2, self.modelljp1.patch.cm)

    def test_Rm_setter(self):
        out = ModelPatch()
        out.Rm = 10e3
        self.assertTrue(np.isclose(out.patch.g_pas, 1/10e3))

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
                         10-self.modelljp2.junction-
                         self.modelljp2.v_low-self.pulse)

    def test_set_vclamp_ls_junction_dur4(self):
        self.assertEqual(self.modelljp2.vclamp.dur4,
                         30)

    def test_set_vclamp_ls_junction_amp5(self):
        self.assertEqual(self.modelljp2.vclamp.amp5,
                         10-self.modelljp2.junction-self.modelljp2.v_low)

    def test_set_vclamp_ls_junction_dur5(self):
        self.assertEqual(self.modelljp2.vclamp.dur5,
                         20)

    def test_set_vclamp_ls_junction_amp6(self):

        self.assertEqual(self.modelljp2.vclamp.amp6,
                         10-self.modelljp2.junction-self.modelljp2.v_low-
                         self.pulse)

    def test_set_vclamp_ls_junction_dur6(self):
        self.assertEqual(self.modelljp2.vclamp.dur6,
                         30)

    def test_set_vclamp_ls_junction_amp7(self):
        self.assertEqual(self.modelljp2.vclamp.amp7,
                         10-self.modelljp2.junction-self.modelljp2.v_low)

    def test_set_vclamp_ls_junction_dur7(self):
        self.assertEqual(self.modelljp2.vclamp.dur7,
                         20)

    def test_set_vclamp_ls_junction_amp8(self):
        self.assertEqual(self.modelljp2.vclamp.amp8,
                         10-self.modelljp2.junction-self.modelljp2.v_low
                         -self.pulse)

    def test_set_vclamp_ls_junction_dur8(self):
        self.assertEqual(self.modelljp2.vclamp.dur8,
                         30)

    def test_set_vclamp_ls_junction_amp9(self):
        self.assertEqual(self.modelljp2.vclamp.amp9,
                         10-self.modelljp2.junction-self.modelljp2.v_low)

    def test_set_vclamp_ls_junction_dur9(self):
        self.assertEqual(self.modelljp2.vclamp.dur9,
                         20)
    def test_set_vclamp_ls_junction_amp10(self):
        self.assertEqual(self.modelljp2.vclamp.amp10,
                         10-self.modelljp2.junction-self.modelljp2.v_low
                         -self.pulse)

    def test_set_vclamp_ls_junction_dur10(self):
        self.assertEqual(self.modelljp2.vclamp.dur10,
                         30)

    def test_set_vclamp_ls_junction_amp11(self):
        self.assertEqual(self.modelljp2.vclamp.amp11,
                         10-self.modelljp2.junction-self.modelljp2.v_low)

    def test_set_vclamp_ls_junction_dur11(self):
        self.assertEqual(self.modelljp2.vclamp.dur11,
                         20)

    def test_set_vclamp_ls_junction_amp12(self):
        self.assertEqual(self.modelljp2.vclamp.amp12,
                         10-self.modelljp2.junction-self.modelljp2.v_low
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



class TestSubtractPassiveProperties(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.modelljp2 = ModelPatch(ljp=10, cm=2)
        cls.dur1 = 10
        cls.dur2 = 20
        cls.delay = 30
        cls.modelljp2.set_vclamp(cls.dur1, 10, cls.dur2, 90, True,
                                 cls.delay)
        cls.pulse = 20
        t_stop = cls.dur1 + 6*cls.dur2 + 5*cls.delay
        cls.current = cls.modelljp2.run(t_stop)


    def test_current_subtraction(self):
        dt = self.modelljp2.dt
        passive = self.current[int(self.dur1/dt)+10:
                               int((self.dur1+self.dur2)/dt)].copy()
        expected = self.modelljp2.curr_stim_response(self.current, self.dur1,
                                                     self.dur2, dt)
        self.assertTrue(np.allclose(passive, expected))

    def test_pulse_sum(self):
        dt = self.modelljp2.dt
        length = int(self.dur2/self.modelljp2.dt)-10
        p_sum = np.zeros((length,))
        t_start = self.dur1 + self.dur2 + 2*self.delay
        for i in range(4):
            p_sum += self.current[int(t_start/dt)+10:
                                  int((t_start+self.dur2)/dt)]
            t_start += self.dur2 + self.delay
        expected = self.modelljp2.curr_leak_amp(self.current, self.dur1,
                                                self.dur2, self.delay, dt)
        self.assertTrue(np.allclose(p_sum, expected))


    def test_leak_substraction(self):
        dt = self.modelljp2.dt
        expected = self.current[int(self.dur1/dt)+10:
                                int((self.dur1+self.dur2)/dt)].copy()
        automatic_leak_subtraction = self.modelljp2.curr_leak_amp(self.current,
                                                                  self.dur1,
                                                                  self.dur2,
                                                                  self.delay,
                                                                  dt)
        
        difference = abs(expected - automatic_leak_subtraction)/expected
        self.assertTrue(np.all(difference < 0.03))

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


        
class TestNonDefaultInit(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.ND = ModelPatch(temp=37, Rm=50000, v_rest=-80,
                            ljp=1, cvode=False,
                            sim_dt=0.01)

    def test_non_default_temp(self):
        self.assertEqual(self.ND.temperature, 37)

    def test_non_default_Rm(self):
        self.assertTrue(np.isclose(self.ND.patch.g_pas, 1/50000))

    def test_non_default_V_rest(self):
        self.assertEqual(self.ND.patch.e_pas, -80)

    def test_non_default_ljp(self):
        self.assertEqual(self.ND.junction, 1)

    def test_non_default_dt(self):
        self.assertEqual(self.ND.sim_dt, 0.01)

    def test_non_default_dt(self):
        self.assertFalse(self.ND.cvode)


if __name__ == "__main__":
    unittest.main()

