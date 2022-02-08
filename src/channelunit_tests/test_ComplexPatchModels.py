import os
import unittest

import numpy as np


from channelunit import ModelWholeCellPatch
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



if __name__ == "__main__":
    unittest.main()
