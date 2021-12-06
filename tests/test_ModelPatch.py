import os
import unittest

import numpy as np

from channelunit import ModelPatch
from channelunit import ModelWholeCellPatch
from channelunit import ModelPatchWithCa


my_loc = os.path.dirname(os.path.abspath(__file__))
channel_loc = os.path.join(my_loc, "..", "demo_CA1", "ion_channels")

F = 96485.33212  # C mol^-1
R = 8.314462618  # J mol^-1 K^-1


class TestModelPatch(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.modelJ = ModelPatch(channel_loc, "nap", "na", 110,
                                gbar_name="gnabar")
        cls.modelNJ = ModelPatch(channel_loc, "nap", "na", 110,
                                 gbar_name="gnabar",
                                 liquid_junction_pot=0)
        cls.modelJ.set_vclamp(10, 10, 100, 100)
        cls.modelNJ.set_vclamp(10, 10, 100, 100)

    def test_reading_in_easy(self):
        out = ModelPatch(channel_loc, "na3", "na", 140)
        self.assertTrue(np.isclose(67.12194015207108, out.E_rev))

    def test_reading_in_no_gbar(self):
        self.assertRaises(SystemExit,  ModelPatch,
                          channel_loc, "na3", "na", gbar_name="gbar1")

    def test_reading_in_provide_E_rev(self):
        out = ModelPatch(channel_loc, "na3", "na", E_rev=40)
        self.assertEqual(40, out.E_rev)

    def test_no_E_rev_name(self):
        self.assertRaises(SystemExit,  ModelPatch,
                          channel_loc, "hd", "nonspecific")

    def test_no_E_rev_name_provide_E_rev(self):
        out = ModelPatch(channel_loc, "hd", "nonspecific", E_rev=-30)
        self.assertEqual(-30, out.E_rev)

    def test_setup_gbar(self):
        out = ModelPatch(channel_loc, "nax", "na", E_rev=-30)
        self.assertEqual(0.001,
                         out.patch.psection()["density_mechs"]["nax"]["gbar"][0])

    def test_setup_ext_conc_e_rev(self):
        out = ModelPatch(channel_loc, "nax", "na", external_conc=140, E_rev=-30)
        conc_fact = np.log(out.external_conc/out.nai)
        new_E_rev = 1e3*R*(273.15+22)/(1*F)*conc_fact
        self.assertEqual(out.E_rev, new_E_rev)



    def test_setup_gbar_custom(self):
        out = ModelPatch(channel_loc, "nap", "na", gbar_name="gnabar")
        self.assertEqual(0.001,
                         out.patch.psection()["density_mechs"]["nap"]["gnabar"][0])

    def test__find_E_rev_value(self):
        out = ModelPatch(channel_loc, "nap", "na", gbar_name="gnabar")
        val = out._find_E_rev_value()
        self.assertEqual(50, val)
        
    
    def test__find_E_rev_name(self):
        out = ModelPatch(channel_loc, "nap", "na", gbar_name="gnabar")
        name = out._find_E_rev_name()
        self.assertEqual("ena", name)

    def test_set_vclamp_junction_amp1(self):
        self.assertEqual(self.modelJ.vclamp.amp1,
                         10-self.modelJ.junction)

    def test_set_vclamp_junction_amp2(self):
        self.assertEqual(self.modelJ.vclamp.amp2,
                         100-self.modelJ.junction)

    def test_set_vclamp_nojunction_amp1(self):
        self.assertEqual(self.modelNJ.vclamp.amp1,
                         10)

    def test_set_vclamp_nojunction_amp2(self):
        self.assertEqual(self.modelNJ.vclamp.amp2,
                         100)

    def test_set_vclamp_junction_dur1(self):
        self.assertEqual(self.modelJ.vclamp.dur1,
                         10)

    def test_set_vclamp_junction_dur2(self):
        self.assertEqual(self.modelJ.vclamp.dur2,
                         100)

    def test_set_vclamp_nojunction_dur1(self):
        self.assertEqual(self.modelNJ.vclamp.dur1,
                         10)

    def test_set_vclamp_nojunction_dur2(self):
        self.assertEqual(self.modelNJ.vclamp.dur2,
                         100)

    def test_extract_current(self):
        I = 2
        out = self.modelJ.extract_current(I, True)
        expected = I/(self.modelJ.vclamp.amp2 - self.modelJ.E_rev)
        self.assertEqual(out, expected)

    def test_extract_current(self):
        I = -22
        out = self.modelJ.extract_current(I, False)
        self.assertEqual(out, 22)

    def test_normalize_to_one(self):
        dic = {1:1, 2:2}
        out = self.modelJ.normalize_to_one(dic)
        self.assertEqual({1:0.5, 2:1.}, out)

    def test_change_nai(self):
        self.modelJ.nai = 5

        conc_fact = np.log(self.modelJ.external_conc/self.modelJ.nai)
        new_E_rev = 1e3*R*(273.15+22)/(1*F)*conc_fact
        self.assertEqual(self.modelJ.E_rev, new_E_rev)

    def test_change_external_conc(self):
        self.modelNJ.external_conc = 140
        conc_fact = np.log(self.modelNJ.external_conc/self.modelNJ.nai)
        new_E_rev = 1e3*R*(273.15+22)/(1*F)*conc_fact
        self.assertEqual(self.modelNJ.E_rev, new_E_rev)
        
class TestCapabilites(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.modelY = ModelPatch(channel_loc, "na3", "na",
                                gbar_name="gbar", E_rev=40, cvode=True)
        cls.stim_levels_act = [-50, -40, -30, -20, -10, 0]
        cls.stim_levels_inact = [-105, -95, -85, -75, -65, -55, -45]
        cls.activationY_cc = cls.modelY.get_activation_steady_state(
            cls.stim_levels_act, -90, 800, 1, chord_conductance=True)
        cls.activationY = cls.modelY.get_activation_steady_state(
            cls.stim_levels_act, -90, 800, 1, chord_conductance=False)
        cls.inactivationY_cc = cls.modelY.get_inactivation_steady_state(
            cls.stim_levels_inact, -5, 20, 1, chord_conductance=True)
        cls.inactivationY = cls.modelY.get_inactivation_steady_state(
            cls.stim_levels_inact, -5, 20, 1, chord_conductance=False)

        cls.modelN = ModelPatch(channel_loc, "na3", "na",
                                gbar_name="gbar", E_rev=40, cvode=False)
        cls.stim_levels_act = [-50, -40, -30, -20, -10, 0]
        cls.stim_levels_inact = [-105, -95, -85, -75, -65, -55, -45]
        cls.activationN_cc = cls.modelN.get_activation_steady_state(
            cls.stim_levels_act, -90, 800, 1, chord_conductance=True)
        cls.activationN = cls.modelN.get_activation_steady_state(
            cls.stim_levels_act, -90, 800, 1, chord_conductance=False)
        cls.inactivationN_cc = cls.modelN.get_inactivation_steady_state(
            cls.stim_levels_inact, -5, 20, 1, chord_conductance=True)
        cls.inactivationN = cls.modelN.get_inactivation_steady_state(
            cls.stim_levels_inact, -5, 20, 1, chord_conductance=False)

    def test_keys_activationY(self):
        self.assertEqual(self.stim_levels_act,
                         list(self.activationY.keys()))

    def test_keys_activationY_cc(self):
        self.assertEqual(self.stim_levels_act,
                         list(self.activationY_cc.keys()))

    def test_keys_inactivationY(self):
        self.assertEqual(self.stim_levels_inact,
                         list(self.inactivationY.keys()))

    def test_keys_inactivationY_cc(self):
        self.assertEqual(self.stim_levels_inact,
                         list(self.inactivationY_cc.keys()))

    def test_keys_activationN(self):
        self.assertEqual(self.stim_levels_act,
                         list(self.activationN.keys()))

    def test_keys_activationN_cc(self):
        self.assertEqual(self.stim_levels_act,
                         list(self.activationN_cc.keys()))

    def test_keys_inactivationN(self):
        self.assertEqual(self.stim_levels_inact,
                         list(self.inactivationN.keys()))

    def test_keys_inactivationN_cc(self):
        self.assertEqual(self.stim_levels_inact,
                         list(self.inactivationN_cc.keys()))
        

class TestWholeCellPatch(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.model = ModelWholeCellPatch(channel_loc, "nap", "na", 110,
                                        gbar_name="gnabar")
        cls.model_init = ModelWholeCellPatch(channel_loc, "nap", "na",
                                             external_conc=110,
                                             gbar_name="gnabar",
                                             temp=35, recompile=False,
                                             liquid_junction_pot=4,
                                             cvode=False,
                                             v_rest=-90,
                                             E_rev=15)
    def test_get_L(self):
        L = self.model.L
        self.assertEqual(L, self.model.patch.L)

    def test_set_L(self):
        self.model.L = 100
        self.assertEqual(100, self.model.patch.L)

    def test_init_ext_conc(self):
        self.assertEqual(self.model_init.external_conc, 110)

    def test_init_temp(self):
        self.assertEqual(self.model_init.temperature, 35)

    def test_E_rev(self):
        conc_fact = np.log(self.model_init.external_conc/self.model_init.nai)
        new_E_rev = 1e3*R*(273.15+self.model_init.temperature)/(1*F)*conc_fact
        self.assertEqual(new_E_rev, self.model_init.E_rev)

    def test_init_ljp(self):
        self.assertEqual(self.model_init.junction, 4)

    def test_init_cvode(self):
        self.assertEqual(self.model_init.cvode, False)

    def test_init_v_rest(self):
        self.assertEqual(self.model_init.patch.e_pas, -90)


class TestPatchWithCa(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.modelcaghk = ModelPatchWithCa(channel_loc, "calHGHK","ca",
                                             external_conc=1.5)
        cls.modelCaghk = ModelPatchWithCa(channel_loc, "CalHGHK","Ca",
                                             external_conc=1.5)
        cls.modelca_eca = ModelPatchWithCa(channel_loc, "calH_eca","ca",
                                              external_conc=1.5)
        cls.modelCa_eCa = ModelPatchWithCa(channel_loc, "CalH_eCa","Ca",
                                              external_conc=1.5)

    def test_raises(self):
        self.assertRaises(SystemExit, ModelPatchWithCa, channel_loc, "callHGHK","cal",
                                             external_conc=1.5)
        
    
if __name__ == "__main__":
    unittest.main()
