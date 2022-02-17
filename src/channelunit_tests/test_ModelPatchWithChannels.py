import os
import unittest

import numpy as np

from channelunit.model_patch import ModelPatchWithChannels
from channelunit import data_path

loc = os.path.dirname(os.path.abspath(__file__))
mechanisms_path = os.path.join(loc, 'mechanisms')
channel_loc = os.path.join(data_path, "ion_channels")
N = 4
DT = 0.1
F = 96485.33212  # C mol^-1
R = 8.314462618  # J mol^-1 K^-1

class TestModelPatchWithChannels(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.modelJ = ModelPatchWithChannels(channel_loc, ["nap"], ["na"],
                                            {"na": 110},
                                            gbar_names={"nap": "gnabar"},
                                            ljp=10, cm=3)

        cls.modelNJ = ModelPatchWithChannels(channel_loc, ["nap"], ["na"],
                                             {"na": 110},
                                             gbar_names={"nap": "gnabar"},
                                             ljp=0)

        cls.modelK = ModelPatchWithChannels(channel_loc, ["kad"], ["k"],
                                            {"k": 2.5},
                                            ljp=0)

        cls.modelca = ModelPatchWithChannels(channel_loc, ["CalH_eCa"], ["Ca"],
                                             {"Ca": 2},
                                             gbar_names={"CalH_eCa":
                                                         "gCalbar"},
                                             gbar_values={"CalH_eCa":0.001},
                                             ljp=0)
        cls.modelmoreions = ModelPatchWithChannels(channel_loc,
                                                   ["nax", "kap"],
                                                   ["na", "k"],
                                                   E_rev={"na": 40})
    def test_cm(self):
        cm = self.modelJ.patch.cm
        self.assertEqual(cm, 3)

    def test_more_than_one_channel(self):
        self.assertEqual(sorted(["nax", "kap"]),
                         sorted(self.modelmoreions.channel_names))

    def test_more_than_one_channel_2(self):
        self.assertEqual(sorted(["na", "k"]),
                         sorted(self.modelmoreions.ion_names))

    def test_more_than_one_channel_3(self):
        self.assertEqual({"na": 40, "k": -77},
                         self.modelmoreions.E_rev)
    def test_reading_in_no_gbar(self):
        self.assertRaises(SystemExit,  ModelPatchWithChannels,
                          channel_loc, ["na3"], ["na"],
                          gbar_names={"na3": "gbar1"})

    def test_setup_gbar(self):
        out = ModelPatchWithChannels(channel_loc, ["nax"], ["na"])
        gbar = out.patch.psection()["density_mechs"]["nax"]["gbar"]
        self.assertEqual(0.001, gbar[0])

    def test_setup_gbar_custom(self):
        out = ModelPatchWithChannels(channel_loc, ["nap"], ["na"],
                                     gbar_names={"nap": "gnabar"})
        gbar = out.patch.psection()["density_mechs"]["nap"]["gnabar"]
        self.assertEqual(0.001, gbar[0])

    def test_setup_gbar_custom_value(self):
        out = ModelPatchWithChannels(channel_loc, ["nap"], ["na"],
                                     gbar_names={"nap": "gnabar"},
                                     gbar_values={"nap": 1})
        gbar = out.patch.psection()["density_mechs"]["nap"]["gnabar"]
        self.assertEqual(1, gbar[0])
        
    def test_Rm(self):
        out = ModelPatchWithChannels(channel_loc, ["nap"], ["na"],
                                     external_conc={"na": 110},
                                     gbar_names={"nap": "gnabar"},
                                     Rm=50000)
        self.assertEqual(out.patch.g_pas, 1/50000)

    def test_gbar_1(self):
        out = ModelPatchWithChannels(channel_loc, ["nap"], ["na"],
                                     external_conc= {"na": 110},
                                     gbar_names={"nap": "gnabar"},
                                     Rm=50000,
                                     gbar_values={"nap": 1})
        self.assertEqual(1, out.get_gbar("nap"))

    def test_gbar_2(self):
        out = ModelPatchWithChannels(channel_loc, ["nap"], ["na"],
                                     external_conc={"na": 110},
                                     gbar_names={"nap": "gnabar"},
                                     Rm=50000, gbar_values={"nap": 1})
        out.set_gbar("nap", 0.2)
        self.assertEqual(0.2, out.get_gbar("nap"))
 
    def test_max_of_dict(self):
        dict1 = {1: np.array([1,3,4,1]), 2: np.array([2,2,2,1])}
        self.assertEqual(ModelPatchWithChannels.get_max_of_dict(dict1),
                         {1: 4, 2:2})

    def test_max_of_dict_1(self):
        dict1 = {1: np.array([-1, -3, -4, -1]), 2: np.array([-2,-2,-2,-1])}
        self.assertEqual(ModelPatchWithChannels.get_max_of_dict(dict1),
                         {1: -4, 2:-2})
    
    def test_max_of_dict_3(self):
        dict1 = {1: np.array([1,3,4,1]), 2: np.array([-2,-2,-2,1])}
        self.assertEqual(ModelPatchWithChannels.get_max_of_dict(dict1),
                         {1: 4, 2:-2})
  
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
        conc_fact = np.log(self.modelJ.external_conc["na"]/5)
        new_E_rev = 1e3*R*(273.15+22)/(1*F)*conc_fact
        self.assertEqual(self.modelJ.E_rev["na"], new_E_rev)

    def test_change_nai_2(self):
        self.modelJ.nai = 5
        conc_fact = np.log(self.modelJ.external_conc["na"]/5)
        new_E_rev = 1e3*R*(273.15+22)/(1*F)*conc_fact
        self.assertEqual(self.modelJ.patch.ena, new_E_rev)

    def test_change_nao(self):
        self.modelNJ.set_external_conc("na", 100)
        conc_fact = np.log(100/self.modelNJ.nai)
        new_E_rev = 1e3*R*(273.15+22)/(1*F)*conc_fact
        self.assertEqual(self.modelNJ.E_rev["na"], new_E_rev)

    def test_change_nai_2(self):
        self.modelNJ.set_external_conc("na", 120)
        conc_fact = np.log(120/self.modelNJ.nai)
        new_E_rev = 1e3*R*(273.15+22)/(1*F)*conc_fact
        self.assertEqual(self.modelNJ.patch.ena, new_E_rev)

    def test_reading_in_easy(self):
        out = ModelPatchWithChannels(channel_loc, ["na3"], ["na"], {"na": 140})
        self.assertTrue(np.isclose(67.12194015207108, out.E_rev["na"]))

    def test_reading_in_easy_2(self):
        out = ModelPatchWithChannels(channel_loc, ["na3"], ["na"], {"na": 140})
        out.run(1)
        self.assertTrue(np.isclose(67.12194015207108, out.patch.ena))
        
    def test_reading_in_provide_E_rev(self):
        out = ModelPatchWithChannels(channel_loc, ["na3"], ["na"], E_rev={"na": 40})
        out.run(1)
        self.assertEqual(40, out.patch.ena)

    def test_reading_in_provide_E_rev_2(self):
        out = ModelPatchWithChannels(channel_loc, ["na3"], ["na"],
                                     E_rev={"na": 40})
        self.assertEqual(40, out.E_rev["na"])

    def test_no_E_rev_name(self):
        self.assertRaises(SystemExit,  ModelPatchWithChannels,
                          channel_loc, ["hd"], ["nonspecific"])

    def test_no_E_rev_name_provide_E_rev(self):
        out = ModelPatchWithChannels(channel_loc, ["hd"], ["nonspecific"],
                                     E_rev={"nonspecific": -30})
        self.assertEqual(-30, out.E_rev["nonspecific"])

    def test_setup_ext_conc_e_rev(self):
        out = ModelPatchWithChannels(channel_loc, ["nax"], ["na"],
                                     external_conc={"na": 140},
                                     E_rev={"na": -30})
        conc_fact = np.log(out.external_conc["na"]/out.nai)
        new_E_rev = 1e3*R*(273.15+22)/(1*F)*conc_fact
        self.assertEqual(out.E_rev["na"], new_E_rev)

    def test_setup_ext_conc_e_rev_2(self):
        out = ModelPatchWithChannels(channel_loc, ["nax"], ["na"],
                                     external_conc={"na": 140},
                                     E_rev={"na": -30})
        conc_fact = np.log(out.external_conc["na"]/out.nai)
        new_E_rev = 1e3*R*(273.15+22)/(1*F)*conc_fact
        out.run(1)
        self.assertEqual(out.patch.ena, new_E_rev)

    def test_calc_E_rev(self):
        out = ModelPatchWithChannels(channel_loc, ["nap"], ["na"],
                                     gbar_names={"nap": "gnabar"})
        val = out.calc_E_rev("na")
        self.assertEqual(50, val) 


class TestAddSingleChannel(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.out = out = ModelPatchWithChannels(channel_loc, ["nax"], ["na"],
                                            {"na": 110})
        cls.out.add_channel("nap", "gnabar", .02)
        
    def test_add_channel_1(self):
        self.assertEqual(["nax", "nap"], self.out.channel_names)

    def test_add_channel_2(self):
        self.assertEqual("gnabar", self.out.gbar_names["nap"])

    def test_add_channel_3(self):
        self.assertEqual([0.02],
                        self.out.patch.psection()["density_mechs"]["nap"]["gnabar"])


class TestCapabilites(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.modelY = ModelPatchWithChannels(channel_loc, ["na3"], ["na"],
                                            gbar_names={"na3": "gbar"},
                                            E_rev={"na": 40},
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

        cls.modelN = ModelPatchWithChannels(channel_loc, ["na3"], ["na"],
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



if __name__ == "__main__":
    unittest.main()
