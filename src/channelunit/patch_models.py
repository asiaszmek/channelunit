from .base_classes import ModelPatch
from .base_classes import ModelPatchCa
from .base_classes import WholeCellAttributes

from .base_classes import memb_shell_width


class ModelWholeCellPatch(ModelPatch, WholeCellAttributes):
    """
    Rin -- in ohms
    """
    def __init__(self, path_to_mods, channel_names: list, ion_names: list,
                 external_conc={}, gbar_names={},
                 temp=22, recompile=True, L=10, diam=10, Ra=100,
                 ljp=0, cvode=True,  Rin=200e6, cap=None, cm=1,
                 v_rest=-65, E_rev={}, gbar_values={}):

        super(ModelWholeCellPatch, self).__init__(path_to_mods,
                                                  channel_names,
                                                  ion_names,
                                                  external_conc=external_conc,
                                                  gbar_names=gbar_names,
                                                  temp=temp,
                                                  recompile=recompile,
                                                  ljp=ljp,
                                                  cvode=cvode, Rm=20000,
                                                  cm=cm, v_rest=v_rest,
                                                  E_rev=E_rev,
                                                  gbar_values=gbar_values)
        self._L = L
        self._diam = diam
        self.patch.L = self._L
        self.patch.Ra = Ra
        self.patch.diam = self._diam
        self._set_g_pas(Rin, [self.patch])
        if cap is None:
            area = self.area([self.patch])*10
            cap = cm*area #should be in nF
            self._cap = cap
        else:
            self.set_cap(cap, [self.patch])



class ModelWholeCellPatchSingleChan(ModelWholeCellPatch):
    def __init__(self, path_to_mods: str, channel_name: str,
                 ion_name: str, external_conc=None,
                 gbar_name="gbar", temp=22, recompile=True,
                 L=10, diam=10, Ra=100,
                 ljp=0, cvode=True,  Rin=200e6,
                 cap=None, cm=1,
                 v_rest=-65, E_rev=None, gbar_value=0.001):
        channel_names = [channel_name]
        ion_names = [ion_name]
        if external_conc is not None:
            ext_conc_dict = {ion_name: external_conc}
        else:
            ext_conc_dict = {}
        gbar_names = {channel_name: gbar_name}
        gbar_values = {channel_name: gbar_value}
        E_rev = {channel_name: E_rev}
        super(ModelWholeCellPatchSingleChan, self).__init__(path_to_mods,
                                                            channel_names,
                                                            ion_names,
                                                            external_conc=ext_conc_dict,
                                                            gbar_names=gbar_names,
                                                            gbar_values=gbar_values,
                                                            E_rev=E_rev,
                                                            temp=temp,
                                                            recompile=recompile,
                                                            ljp=ljp,
                                                            cvode=cvode,
                                                            v_rest=v_rest,
                                                            Rin=Rin, cap=cap,
                                                            cm=cm, L=L,
                                                            diam=diam, Ra=Ra)
        


class ModelWholeCellPatchCa(ModelPatchCa, WholeCellAttributes):
    def __init__(self, path_to_mods, channel_names: list, ion_names: list,
                 external_conc: dict, E_rev={}, gbar_names={}, temp=22, recompile=True,
                 ljp=0, cvode=True, Rin=200e6, cap=None, cm=1,
                 v_rest=-65,
                 gbar_values={}, t_decay=20, L=10, diam=10, Ra=100,
                 buffer_capacity=18,
                 membrane_shell_width=memb_shell_width):
        """
        Model class for testing calcium channels with 
        cap -- nF, cell capacitance 
        cm -- membrane capacitance 1 uF/cm2
        """
            
        super(ModelWholeCellPatchCa, self).__init__(path_to_mods,
                                                    channel_names,
                                                    ion_names,
                                                    external_conc=external_conc,
                                                    E_rev=E_rev,
                                                    gbar_names=gbar_names,
                                                    temp=temp,
                                                    recompile=recompile,
                                                    ljp=ljp,
                                                    cvode=cvode,
                                                    Rm=20000, cm=cm,
                                                    v_rest=v_rest,
                                                    gbar_values=gbar_values,
                                                    t_decay=t_decay, Ra=Ra,
                                                    buffer_capacity=buffer_capacity,
                                                    membrane_shell_width=membrane_shell_width)
        self._L = L
        self._diam = diam
        self.patch.L = self._L
        self.patch.Ra = Ra
        self.patch.diam = self._diam
        self._set_g_pas(Rin, [self.patch])
        if cap is None:
            area = self.area([self.patch])*10
            cap = cm*area #should be in nF
            self._cap = cap
        else:
            self.set_cap(cap, [self.patch])


class ModelWholeCellPatchCaSingleChan(ModelWholeCellPatchCa):
    def __init__(self, path_to_mods: str, channel_name: str, ion_name: str,
                 external_conc,
                 E_rev={}, gbar_name="gbar", temp=22, recompile=True, L=10,
                 diam=10, Ra=100,
                 ljp=0, cvode=True,  Rin=200e6, cap=None, cm=1,
                 v_rest=-65, gbar_value=0.001,
                 t_decay=20,
                 buffer_capacity=18,
                 membrane_shell_width=memb_shell_width):
        channel_names = [channel_name]
        ion_names = [ion_name]
        if external_conc is not None:
            ext_conc_dict = {ion_name: external_conc}
        else:
            raise SystemError("Ca external conc needs to be specified")
        gbar_names = {channel_name: gbar_name}
        gbar_values = {channel_name: gbar_value}
        super(ModelWholeCellPatchCaSingleChan, self).__init__(path_to_mods,
                                                              channel_names,
                                                              ion_names,
                                                              external_conc=ext_conc_dict,
                                                              E_rev=E_rev,
                                                              gbar_names=gbar_names,
                                                              gbar_values=gbar_values,
                                                              temp=temp,
                                                              recompile=recompile,
                                                              ljp=ljp,
                                                              cvode=cvode,
                                                              v_rest=v_rest,
                                                              Rin=Rin,
                                                              cap=cap, cm=cm,
                                                              t_decay=t_decay,
                                                              L=L, diam=diam,
                                                              Ra=Ra,
                                                              buffer_capacity=buffer_capacity,
                                                              membrane_shell_width=membrane_shell_width)


class ModelOocytePatch(ModelWholeCellPatch):
    #parameters from  PMID: 20737886 DOI: 10.1016/0012-1606(81)90417-6 
    def __init__(self, path_to_mods, channel_names: list, ion_names: list,
                 external_conc={}, gbar_names={},
                 temp=22, recompile=True, ljp=0, cvode=True,
                 E_rev={}, gbar_values={}):
        super(ModelOocyte, self).__init__(path_to_mods, channel_names,
                                          ion_names,
                                          external_conc=external_conc,
                                          gbar_names=gbar_names,
                                          temp=22, recompile=True, L=1.3e3,
                                          diam=1.3e3, Ra=100,
                                          ljp=0, cvode=True,
                                          Rin=1.86e6, cap=None, cm=12,
                                          v_rest=-50, E_rev=E_rev,
                                          gbar_values=gbar_values)


class ModelOocytePatchCa(ModelWholeCellPatchCa):
    #parameters from  PMID: 20737886 DOI: 10.1016/0012-1606(81)90417-6
    # (electric)
    #Ca params from https://doi.org/10.1016/j.ydbio.2005.10.034
    def __init__(self, path_to_mods, channel_names: list, ion_names: list,
                 external_conc={}, E_rev={}, gbar_names={},
                 temp=22, recompile=True, ljp=0, cvode=True,
                 gbar_values={}):
        super(ModelOocyteCa, self).__init__(path_to_mods, channel_names,
                                            ion_names,
                                            external_conc=external_conc,
                                            E_rev=E_rev, gbar_names=gbar_names,
                                            temp=22, recompile=True,
                                            L=1.3e3, diam=1.3e3, Ra=100,
                                            ljp=0, cvode=True,
                                            Rin=1.86e6, cap=None, cm=12,
                                            v_rest=-50, gbar_values=gbar_values,
                                            t_decay=8e3)


class ModelGiantExcisedPatch(ModelPatch):
    pass


class ModelGiantExcisedPatchCa(ModelPatchCa):
    pass
