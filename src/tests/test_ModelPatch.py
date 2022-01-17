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
        cls.modelljp2.set_vclamp(10, 10, 20, 90, True, 30)
        cls.pulse = 20

    def test_set_vclamp_junction_amp1(self):
        self.assertEqual(self.modelljp1.vclamp.amp1,
                         10-self.modelljp1.junction)

    def test_set_vclamp_junction_amp2(self):
        self.assertEqual(self.modelljp1.vclamp.amp2,
                         100-self.modelljp1.junction)
                                                                                                                                                                          
    def test_set_vclamp_junction_dur1(self):                                                                                                                              
        self.assertEqual(self.modelljp1.vclamp.dur1,                                                                                                                      
                         10)                                                                                                                                              
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
        print(self.modelljp2.vclamp.amp4)
        self.assertEqual(self.modelljp2.vclamp.amp4,
                         10-2*self.modelljp2.junction-self.modelljp2.v_low-self.pulse)

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
                         10-2*self.modelljp2.junction-self.modelljp2.v_low-self.pulse)

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
                         self.modelljp2.vclamp.amp4+self.pulse/4)

    def test_pulse_height_2(self):
        self.assertEqual(self.modelljp2.vclamp.amp7,
                         self.modelljp2.vclamp.amp6+self.pulse/4)

    def test_pulse_height_3(self):                                                                                                                                        
        self.assertEqual(self.modelljp2.vclamp.amp9,                                                                                                                      
                         self.modelljp2.vclamp.amp8+self.pulse/4)                                                                                                         
                                                                                                                                                                          
    def test_pulse_height_4(self):                                                                                                                                        
        self.assertEqual(self.modelljp2.vclamp.amp11,                                                                                                                     
                         self.modelljp2.vclamp.amp10+self.pulse/4)                                                                                                        
                                                                                                                                                                          
    def test_leak_subtraction(self):
        pass


if __name__ == "__main__":
    unittest.main()

