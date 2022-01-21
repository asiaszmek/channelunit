import os
import unittest

import numpy as np


from channelunit.model_patch import ModelPatchWithChannels
from channelunit import ModelWholeCellPatch
from channelunit import ModelWholeCellPatchCaShell
from channelunit import ModelWholeCellPatchCaShellOneChannel
from channelunit import data_path



channel_loc = os.path.join(data_path, "ion_channels")

F = 96485.33212  # C mol^-1
R = 8.314462618  # J mol^-1 K^-1

class TestModelPatchWithChannels(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.modelJ = ModelPatchWithChannels(channel_loc, ["nap"], ["na"], {"na": 110},
                                gbar_names={"nap": "gnabar"},
                                liquid_junction_pot=10)
        
       
        cls.modelJ2 = ModelPatchWithChannels(channel_loc, ["nap"], ["na"], {"na": 110},
                                 gbar_names={"nap": "gnabar"},
                                 liquid_junction_pot=10)

    def test_reading_in_no_gbar(self):
        self.assertRaises(SystemExit,  ModelPatchWithChannels,
                          channel_loc, ["na3"], ["na"], gbar_names={"na3": "gbar1"})

    def test_setup_gbar(self):
        out = ModelPatchWithChannels(channel_loc, ["nax"], ["na"])
        self.assertEqual(0.001,
                         out.patch.psection()["density_mechs"]["nax"]["gbar"][0])

    def test_setup_gbar_custom(self):
        out = ModelPatchWithChannels(channel_loc, ["nap"], ["na"], gbar_names={"nap": "gnabar"})
        self.assertEqual(0.001,
                         out.patch.psection()["density_mechs"]["nap"]["gnabar"][0])

    def test_setup_gbar_custom_value(self):
        out = ModelPatchWithChannels(channel_loc, ["nap"], ["na"], gbar_names={"nap": "gnabar"},
                                     gbar_values={"nap": 1})
        self.assertEqual(1,
                         out.patch.psection()["density_mechs"]["nap"]["gnabar"][0])
        
    def test_Rm(self):
        out = ModelPatchWithChannels(channel_loc, ["nap"], ["na"],
                                     external_conc={"na": 110},
                                     gbar_names={"nap": "gnabar"},
                                     R_m=50000)
        self.assertEqual(out.patch.g_pas, 1/50000)

    def test_gbar_1(self):
        out = ModelPatchWithChannels(channel_loc, ["nap"], ["na"],
                                     external_conc= {"na": 110},
                                     gbar_names={"nap": "gnabar"},
                                     R_m=50000,
                                     gbar_values={"nap": 1})
        self.assertEqual(1, out.get_gbar("nap"))

    def test_gbar_2(self):
        out = ModelPatchWithChannels(channel_loc, ["nap"], ["na"],
                                     external_conc={"na": 110},
                                     gbar_names={"nap": "gnabar"},
                                     R_m=50000, gbar_values={"nap": 1})
        out.set_gbar("nap", 0.2)
        self.assertEqual(0.2, out.get_gbar("nap"))
 
    def test_max_of_dict(self):
        dict1 = {1: np.array([1,3,4,1]), 2: np.array([2,2,2,1])}
        self.assertEqual(ModelPatchWithChannels.get_max_of_dict(dict1, False, ["k"], False),
                         {1: 4, 2:2})

    def test_max_of_dict_1(self):
        dict1 = {1: np.array([-1, -3, -4, -1]), 2: np.array([-2,-2,-2,-1])}
        self.assertEqual(ModelPatchWithChannels.get_max_of_dict(dict1, False, ["na"], False),
                         {1: -4, 2:-2})
    
    def test_max_of_dict_3(self):
        dict1 = {1: np.array([1,3,4,1]), 2: np.array([2,2,2,1])}
        self.assertEqual(ModelPatchWithChannels.get_max_of_dict(dict1, True, ["na"], True),
                         {1: 4, 2:2})

        
class TestModelWholeCell(unittest.TestCase):
    @classmethod
    def setUpClass(cls):

        cls.modelJ = ModelWholeCellPatch(channel_loc, ["nap"], ["na"], {"na": 110},
                                         gbar_names={"nap": "gnabar"})
        cls.modelNJ = ModelWholeCellPatch(channel_loc, ["nap"], ["na"], {"na": 110},
                                          gbar_names={"nap": "gnabar"},
                                          liquid_junction_pot=0)
    
    def test_reading_in_easy(self):
        out = ModelWholeCellPatch(channel_loc, ["na3"], ["na"],{"na": 140})
        self.assertTrue(np.isclose(67.12194015207108, out.E_rev["na"]))

    def test_reading_in_provide_E_rev(self):
        out = ModelWholeCellPatch(channel_loc, ["na3"], ["na"], E_rev={"na": 40})
        self.assertEqual(40, out.E_rev["na"])

    def test_no_E_rev_name(self):
        self.assertRaises(SystemExit,  ModelWholeCellPatch,
                          channel_loc, ["hd"], ["nonspecific"])

    def test_no_E_rev_name_provide_E_rev(self):
        out = ModelWholeCellPatch(channel_loc, ["hd"], ["nonspecific"],
                                  E_rev={"nonspecific": -30})
        self.assertEqual(-30, out.E_rev["nonspecific"])

    def test_setup_ext_conc_e_rev(self):
        out = ModelWholeCellPatch(channel_loc, ["nax"], ["na"],
                                  external_conc={"na": 140},
                                  E_rev={"na": -30})
        conc_fact = np.log(out.external_conc["na"]/out.nai)
        new_E_rev = 1e3*R*(273.15+22)/(1*F)*conc_fact
        self.assertEqual(out.E_rev["na"], new_E_rev)

    def testcalc_E_rev(self):
        out = ModelWholeCellPatch(channel_loc, ["nap"], ["na"], gbar_names={"nap": "gnabar"})
        val = out.calc_E_rev("na")
        self.assertEqual(50, val) 

    def test_extract_current_chord_conductance(self):
        dt = 0.01
        I = np.ones((int((10+20)/dt)))
        out = self.modelJ.extract_current(I, True, False, 10, 20, 0, dt)
        expected = I[int(10/dt)+10:]/(self.modelJ.vclamp.amp2
                                       - self.modelJ.E_rev["na"])
        comparison = np.allclose(expected, out)
        self.assertTrue(comparison)

    def test_extract_current_nothing(self):
        dt = 0.01
        I = np.ones((int((10+20)/dt)))
        out = self.modelJ.extract_current(I, False, False, 10, 20, 0, dt)
        expected = I[int(10/dt)+10:]
        comparison = np.allclose(expected, out)
        self.assertTrue(comparison)

    def test_extract_current_ls_chord_conductance(self):
        dt = 0.01
        dur1 = 20
        dur2 = 10
        delay = 5
        I = np.ones((int((dur1+5*dur2+6*delay)/dt)))
        I[:int(dur1/dt)+10] = 0
        I[int(dur1/dt):int((dur1+dur2)/dt)] = 10
        I[int((dur1+dur2)/dt):int((dur1+dur2+delay)/dt)] = 0
        I[int((dur1+dur2+delay)/dt):int((dur1+dur2+2*delay)/dt)] = -5
        t_start = dur1+dur2+2*delay
        for i in range(4):
            I[int(t_start/dt):int((t_start+dur2)/dt)] = 2
            I[int((t_start+dur2)/dt):int((t_start+dur2+delay)/dt)] = -5
            t_start += dur2 + delay
        out = self.modelJ.extract_current(I, True, True, dur1, dur2, delay,
                                          dt)
        expected = (I[int(dur1/dt)+10:int((dur1+dur2)/dt)]
                    -8)/(self.modelJ.vclamp.amp2 - self.modelJ.E_rev["na"])
        comparison = np.allclose(expected, out)
        self.assertTrue(comparison)

    def test_extract_current_ls_no_chord_conductance(self):
        dt = 0.01
        dur1 = 20
        dur2 = 10
        delay = 5
        I = np.ones((int((dur1+5*dur2+6*delay)/dt)))
        I[:int(dur1/dt)] = 0
        I[int(dur1/dt):int((dur1+dur2)/dt)] = 10
        I[int((dur1+dur2)/dt):int((dur1+dur2+delay)/dt)] = 0
        I[int((dur1+dur2+delay)/dt):int((dur1+dur2+2*delay)/dt)] = -5
        t_start = dur1+dur2+2*delay
        for i in range(4):
            I[int(t_start/dt):int((t_start+dur2)/dt)] = 2
            I[int((t_start+dur2)/dt):int((t_start+dur2+delay)/dt)] = -5
            t_start += dur2 + delay
        out = self.modelJ.extract_current(I, False, True, dur1, dur2, delay,
                                          dt)
        expected = (I[int(dur1/dt)+10:int((dur1+dur2)/dt)]-8)
        comparison = np.allclose(expected, out)
        self.assertTrue(comparison)

        
    def test_normalize_to_one(self):
        dic = {1:-1, 2:-2}
        out = self.modelJ.normalize_to_one(dic, "save_sign")
        self.assertEqual({1:-0.5, 2:-1.}, out)

    def test_normalize_to_one_2(self):
        dic = {1:-1, 2:-2}
        out = self.modelJ.normalize_to_one(dic, "to_one")
        self.assertEqual({1:0.5, 2:1.}, out)

    def test_change_nai(self):
        self.modelJ.nai = 5
        conc_fact = np.log(self.modelJ.external_conc["na"]/self.modelJ.nai)
        new_E_rev = 1e3*R*(273.15+22)/(1*F)*conc_fact
        self.assertEqual(self.modelJ.E_rev["na"], new_E_rev)


    def test_setup_gbar_custom_value(self):
        out = ModelWholeCellPatch(channel_loc, ["nap"], ["na"], gbar_names={"nap": "gnabar"},
                                  gbar_values={"nap": 1}, E_rev={"na": 50})
        self.assertEqual(1,
                         out.patch.psection()["density_mechs"]["nap"]["gnabar"][0])

        
class TestCapabilites(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.modelY = ModelWholeCellPatch(channel_loc, ["na3"], ["na"],
                                         gbar_names={"na3": "gbar"}, E_rev={"na": 40},
                                         cvode=True)
        cls.stim_levels_act = [-50, -40, -30, -20, -10, 0]
        cls.stim_levels_inact = [-105, -95, -85, -75, -65, -55, -45]
        cls.activationY_cc = cls.modelY.get_activation_steady_state(
            cls.stim_levels_act, -90, 200, 1, chord_conductance=True,
            electrode_current=True)
        cls.activationY = cls.modelY.get_activation_steady_state(
            cls.stim_levels_act, -90, 200, 1, chord_conductance=False,
            electrode_current=True)
        cls.inactivationY_cc = cls.modelY.get_inactivation_steady_state(
            cls.stim_levels_inact, -5, 20, 1, chord_conductance=True,
            electrode_current=True)
        cls.inactivationY = cls.modelY.get_inactivation_steady_state(
            cls.stim_levels_inact, -5, 20, 1, chord_conductance=False,
            electrode_current=True)

        cls.modelN = ModelWholeCellPatch(channel_loc, ["na3"], ["na"],
                                         gbar_names={"na3": "gbar"},
                                         E_rev={"na": 40},
                                         cvode=False)
        cls.stim_levels_act = [-50, -40, -30, -20, -10, 0]
        cls.stim_levels_inact = [-105, -95, -85, -75, -65, -55, -45]
        cls.activationN_cc = cls.modelN.get_activation_steady_state(
            cls.stim_levels_act, -90, 200, 1, chord_conductance=True,
            electrode_current=True)
        cls.activationN = cls.modelN.get_activation_steady_state(
            cls.stim_levels_act, -90, 200, 1, chord_conductance=False,
            electrode_current=True)
        cls.inactivationN_cc = cls.modelN.get_inactivation_steady_state(
            cls.stim_levels_inact, -5, 20, 1, chord_conductance=True,
            electrode_current=True)
        cls.inactivationN = cls.modelN.get_inactivation_steady_state(
            cls.stim_levels_inact, -5, 20, 1, chord_conductance=False,
            electrode_current=True)

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
        

class TestWholeCellPatchNernst(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.model = ModelWholeCellPatch(channel_loc, ["nap"], ["na"], {"na": 110},
                                        gbar_names={"nap": "gnabar"})
        cls.model_init = ModelWholeCellPatch(channel_loc, ["nap"], ["na"],
                                             external_conc={"na": 110},
                                             gbar_names={"nap": "gnabar"},
                                             temp=35, recompile=False,
                                             liquid_junction_pot=4,
                                             cvode=False,
                                             v_rest=-90,
                                             E_rev={"na": 15}, R_in=2e9)
    def test_setup_gbar_custom_value(self):
        out = ModelWholeCellPatch(channel_loc, ["nap"], ["na"],
                                  gbar_names={"nap": "gnabar"},
                                  gbar_values={"nap": 1},
                                  external_conc={"na": 110},
                                  temp=35, recompile=False,
                                  liquid_junction_pot=4,
                                  cvode=False,
                                  v_rest=-90,
                                  E_rev={"na": 15}, R_in=2e9)
        self.assertEqual(1,
                         out.patch.psection()["density_mechs"]["nap"]["gnabar"][0])

    def test_get_L(self):
        L = self.model.L
        self.assertEqual(L, self.model.patch.L)

    def test_set_L(self):
        self.model.L = 100
        self.assertEqual(100, self.model.patch.L)

    def test_init_ext_conc(self):
        self.assertEqual(self.model_init.external_conc["na"], 110)

    def test_init_temp(self):
        self.assertEqual(self.model_init.temperature, 35)

    def test_E_rev(self):
        conc_fact = np.log(self.model_init.external_conc["na"]/self.model_init.nai)
        new_E_rev = 1e3*R*(273.15+self.model_init.temperature)/(1*F)*conc_fact
        self.assertEqual(new_E_rev, self.model_init.E_rev["na"])

    def test_init_ljp(self):
        self.assertEqual(self.model_init.junction, 4)

    def test_init_cvode(self):
        self.assertEqual(self.model_init.cvode, False)

    def test_init_v_rest(self):
        self.assertEqual(self.model_init.patch.e_pas, -90)
    
    def test_set_R_in(self):
        area = 0
        for seg in self.model_init.patch:
            area += seg.area()*1e-4
        self.assertEqual(self.model_init.patch.g_pas, 1/(area)/2*1e-9)

class TestPatchWithCa(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.modelcaghk = ModelWholeCellPatchCaShellOneChannel(channel_loc, "calHGHK", "ca",
                                                   1.5)
        cls.modelCaghk = ModelWholeCellPatchCaShellOneChannel(channel_loc, "CalHGHK","Ca",
                                                   1.5)
        cls.modelca_eca = ModelWholeCellPatchCaShellOneChannel(channel_loc, "calH_eca",
                                                    "ca",
                                                    1.5,
                                                    gbar_name="gcal")
        cls.modelCa_eCa = ModelWholeCellPatchCaShellOneChannel(channel_loc, "CalH_eCa",
                                                    "Ca", 1.5, gbar_name="gCal")

    def test_raises(self):
        self.assertRaises(SystemExit, ModelWholeCellPatchCaShellOneChannel, channel_loc,
                          "callHGHK","cal", 1.5)





if __name__ == "__main__":
    unittest.main()
