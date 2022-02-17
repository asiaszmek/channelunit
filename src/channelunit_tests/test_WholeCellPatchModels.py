import os
import unittest

import numpy as np


from channelunit import ModelWholeCellPatch
from channelunit import ModelWholeCellPatchOneChannel
from channelunit import ModelWholeCellPatchCaShell
from channelunit import ModelWholeCellPatchCaShellOneChannel
from channelunit import data_path

channel_loc = os.path.join(data_path, "ion_channels")

F = 96485.33212  # C mol^-1
R = 8.314462618  # J mol^-1 K^-1
       

class TestModelWholeCell(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.model = ModelWholeCellPatch(channel_loc, ["nap"], ["na"],
                                        {"na": 110},
                                        gbar_names={"nap": "gnabar"}, cap=1, cm=2)
        cls.model_init = ModelWholeCellPatch(channel_loc, ["nap"], ["na"],
                                             external_conc={"na": 110},
                                             gbar_names={"nap": "gnabar"},
                                             temp=35, recompile=False,
                                             ljp=4,
                                             cvode=False,
                                             v_rest=-90,
                                             E_rev={"na": 15}, Rin=2e9, cm=2)

    def test_cap_init(self):
        area = self.model.area([self.model.patch])*1e-8
        expected = 1/area*1e5
        self.assertEqual(self.model.cm, expected)

    def test_cm_init(self):
        self.assertEqual(2, self.model_init.patch.cm)
        
    def test_setup_gbar_custom_value(self):
        out = ModelWholeCellPatch(channel_loc, ["nap"], ["na"],
                                  gbar_names={"nap": "gnabar"},
                                  gbar_values={"nap": 1},
                                  external_conc={"na": 110},
                                  temp=35, recompile=False,
                                  ljp=4,
                                  cvode=False,
                                  v_rest=-90,
                                  E_rev={"na": 15}, Rin=2e9)
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
        conc_fact = np.log(self.model_init.external_conc["na"]
                           /self.model_init.nai)
        new_E_rev = 1e3*R*(273.15+self.model_init.temperature)/(1*F)*conc_fact
        self.assertEqual(new_E_rev, self.model_init.E_rev["na"])

    def test_init_ljp(self):
        self.assertEqual(self.model_init.junction, 4)

    def test_init_cvode(self):
        self.assertEqual(self.model_init.cvode, False)

    def test_init_v_rest(self):
        self.assertEqual(self.model_init.patch.e_pas, -90)
    
    def test_set_Rin(self):
        area = 0
        for seg in self.model_init.patch:
            area += seg.area()*1e-8
        self.assertTrue(np.isclose(self.model_init.patch.g_pas, 1/(area)/2*1e-9))

    def test_setup_cap(self):
        out = ModelWholeCellPatch(channel_loc, ["nap"], ["na"],
                                  {"na": 110},
                                  gbar_names={"nap": "gnabar"})
        out.cap = 2
        area = out.area([out.patch])*1e-8
        self.assertEqual(out.patch.cm, 2/area*1e5)
    

class TestChangeL(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.model = ModelWholeCellPatch(channel_loc, ["nap"], ["na"],
                                        {"na": 110},
                                        gbar_names={"nap": "gnabar"}, cap=1, cm=2)
        cls.old_gpas = cls.model.patch.g_pas
        cls.model.L = 100

    def test_changing_L(self):
        self.assertEqual(self.model.L, 100)


    def test_same_diam(self):
        self.assertEqual(self.model.diam, 10)

    def test_changing_g_pas(self):
        self.assertTrue(np.isclose(self.model.patch.g_pas, self.old_gpas/10))


class TestChangeDiam(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.model = ModelWholeCellPatch(channel_loc, ["nap"], ["na"],
                                        {"na": 110},
                                        gbar_names={"nap": "gnabar"}, cap=1, cm=2)
        cls.old_gpas = cls.model.patch.g_pas
        cls.model.diam = 100

    def test_changing_diam(self):
        self.assertEqual(self.model.diam, 100)

    def test_same_L(self):
        self.assertEqual(self.model.L, 10)

    def test_changing_g_pas(self):
        self.assertEqual(self.model.patch.g_pas, self.old_gpas/10)

        
class TestPatchWithCa(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.modelcaghk = ModelWholeCellPatchCaShellOneChannel(channel_loc,
                                                              "calHGHK", "ca",
                                                              1.5)
        cls.modelCaghk = ModelWholeCellPatchCaShellOneChannel(channel_loc,
                                                              "CalHGHK","Ca",
                                                              1.5)
        cls.modelca_eca = ModelWholeCellPatchCaShellOneChannel(channel_loc,
                                                               "calH_eca",
                                                               "ca",
                                                               1.5,
                                                               gbar_name="gcal")
        cls.modelCa_eCa = ModelWholeCellPatchCaShellOneChannel(channel_loc,
                                                               "CalH_eCa",
                                                               "Ca", 1.5,
                                                               gbar_name="gCal")

    def test_raises(self):
        self.assertRaises(SystemExit, ModelWholeCellPatchCaShellOneChannel,
                          channel_loc,
                          "callHGHK", "cal", 1.5)


class TestModelWholeCellPatchOneChannel(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.model = ModelWholeCellPatchOneChannel(channel_loc,
                                                  "nap", "na", external_conc=110,
                                                  gbar_name="gnabar")

    def test_channel_names(self):
        self.assertEqual(self.model.channel_names, ["nap"])

    def test_channels(self):
        keys = self.model.patch.psection()["density_mechs"]
        self.assertEqual(sorted(keys), ["nap", "pas"])

    def test_ion_names(self):
        self.assertEqual(["na"], self.model.ion_names)

    def test_ext_conc(self):
        self.assertEqual(self.model.external_conc, {"na": 110, "Ca": None})

    def test_gbar_name(self):
        self.assertEqual(self.model.gbar_names, {"nap": "gnabar"}) 


if __name__ == "__main__":
    unittest.main()
