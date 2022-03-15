import unittest
import numpy as np

from channelunit.base_classes import MembranePatch


class TestVclamp(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.modelljp1 = MembranePatch(temp=22, Rm=20000, cm=1,
                                      v_rest=-65, ljp=10,
                                      cvode=True, sim_dt=0.001)
        cls.modelljp2 = MembranePatch(temp=22, Rm=20000, cm=2,
                                      v_rest=-65, ljp=10,
                                      cvode=False, sim_dt=0.001)
        cls.modelljp1.set_vclamp(10, 10, 100, 100, False)
        cls.dur1 = 10
        cls.dur2 = 20
        cls.delay = cls.dur1
        cls.modelljp2.set_vclamp(cls.dur1, 10, cls.dur2, 90, True)
        cls.pulse = 20

    def test_sim_dt(self):
        self.assertEqual(self.modelljp2.sim_dt, 0.001)

    def test_init_1(self):
        self.assertEqual(self.modelljp1.temperature, 22)

    def test_init_pas_1(self):
        self.assertTrue(np.isclose(self.modelljp1.patch.g_pas,
                                   1/20000))

    def test_init_pas_2(self):
        self.assertTrue(np.isclose(self.modelljp1.patch.e_pas,
                                   -65))

    def test_init_cvode(self):
        self.assertTrue(self.modelljp1.cvode)

    def test_rs(self):
        self.modelljp1.vclamp.rs = 100
        self.assertEqual(100, self.modelljp1.vclamp.rs)
    
    def test_getter_cm(self):
        cm = self.modelljp2.patch.cm
        self.assertEqual(cm, 2)

    def test_setter_cm(self):
        self.modelljp1.cm = 2
        self.assertEqual(2, self.modelljp1.patch.cm)

    def test_Rm_setter(self):
        self.modelljp2.Rm = 10e3
        self.assertTrue(np.isclose(self.modelljp2.patch.g_pas, 1/10e3))

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
                         20)

    def test_set_vclamp_ls_junction_amp4(self):
        self.assertEqual(self.modelljp2.vclamp.amp4,
                         10-self.modelljp2.junction-
                         self.modelljp2.v_low-self.pulse)

    def test_set_vclamp_ls_junction_dur4(self):
        self.assertEqual(self.modelljp2.vclamp.dur4,
                         20)

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
                         20)

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
                         20)

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
                         20)

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
                         20)

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




if __name__ == "__main__":
    unittest.main()

