from .base_classes import ModelPatch
from .base_classes import ModelPatchCa
from .base_classes import WholeCellAttributes

memb_shell_width = .1

class ModelWholeCellPatch(ModelPatch, WholeCellAttributes):
    """
    Rin -- in ohms
    """
    def __init__(self, path_to_mods, channel_names: list, ion_names: list,
                 external_conc={}, internal_conc={}, gbar_names={},
                 temp=22, recompile=True, L=10, diam=10, Ra=100,
                 ljp=0, cvode=True,  Rin=200e6, Rm=20000, cap=None, cm=1,
                 v_rest=-65, E_rev={}, gbar_values={}):

        super(ModelWholeCellPatch, self).__init__(path_to_mods,
                                                  channel_names,
                                                  ion_names,
                                                  external_conc=external_conc,
                                                  internal_conc=internal_conc,
                                                  gbar_names=gbar_names,
                                                  temp=temp,
                                                  recompile=recompile,
                                                  ljp=ljp,
                                                  cvode=cvode, Rm=Rm,
                                                  cm=cm, v_rest=v_rest,
                                                  E_rev=E_rev,
                                                  gbar_values=gbar_values,
                                                  sim_dt=0.001,
                                                  directory="")
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
        if Rin is None:
            area = self.area([self.patch])*1e-8
            self._Rin = area * Rm
        else:
            self._set_g_pas(Rin, [self.patch])
        self.vclamp.rs = 10

class ModelWholeCellPatchSingleChan(ModelWholeCellPatch):
    def __init__(self, path_to_mods: str, channel_name: str,
                 ion_name: str, external_conc=None,
                 internal_conc=None,
                 gbar_name="gbar", temp=22, recompile=True,
                 L=10, diam=10, Ra=100,
                 ljp=0, cvode=True,  Rin=200e6, Rm=20000,
                 cap=None, cm=1,
                 v_rest=-65, E_rev=None, gbar_value=0.001):
        channel_names = [channel_name]
        ion_names = [ion_name]
        ext_c = {}
        int_c = {}     
        if external_conc is not None:
            ext_c[ion_name] = external_conc
        if internal_conc is not None:
            int_c[ion_name] = internal_conc
        gbar_n = {channel_name: gbar_name}
        gbar_v = {channel_name: gbar_value}
        E_rev = {channel_name: E_rev}
        super(ModelWholeCellPatchSingleChan, self).__init__(path_to_mods,
                                                            channel_names,
                                                            ion_names,
                                                            external_conc=ext_c,
                                                            internal_conc=int_c,
                                                            gbar_names=gbar_n,
                                                            gbar_values=gbar_v,
                                                            E_rev=E_rev,
                                                            temp=temp,
                                                            recompile=recompile,
                                                            ljp=ljp,
                                                            cvode=cvode,
                                                            v_rest=v_rest,
                                                            Rin=Rin, Rm=Rm,
                                                            cap=cap,
                                                            cm=cm, L=L,
                                                            diam=diam, Ra=Ra)
        


class ModelWholeCellPatchCa(ModelPatchCa, WholeCellAttributes):
    def __init__(self, path_to_mods, channel_names: list, ion_names: list,
                 external_conc: dict, internal_conc={}, E_rev={},
                 gbar_names={}, temp=22, recompile=True,
                 ljp=0, cvode=True, Rin=200e6, Rm=20000, cap=None, cm=1,
                 v_rest=-65,
                 gbar_values={}, t_decay=20, L=10, diam=10, Ra=100,
                 buffer_capacity=18,
                 membrane_shell_width=memb_shell_width,
                 directory="validation_results", sim_dt=0.001):
        """
        Model class for testing calcium channels with 
        cap -- nF, cell capacitance 
        cm -- membrane capacitance 1 uF/cm2
        """
            
        super(ModelWholeCellPatchCa, self).__init__(path_to_mods,
                                                    channel_names,
                                                    ion_names,
                                                    external_conc=external_conc,
                                                    internal_conc=internal_conc,
                                                    E_rev=E_rev,
                                                    gbar_names=gbar_names,
                                                    temp=temp,
                                                    recompile=recompile,
                                                    ljp=ljp,
                                                    cvode=cvode, L=L, diam=diam,
                                                    Rm=Rm, cm=cm,
                                                    v_rest=v_rest,
                                                    gbar_values=gbar_values,
                                                    t_decay=t_decay, Ra=Ra,
                                                    buffer_capacity=buffer_capacity,
                                                    membrane_shell_width=membrane_shell_width,
                                                    directory=directory,
                                                    sim_dt=sim_dt)
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
        if Rin is None:
            area = self.area([self.patch])*1e-8
            self._Rin = area * Rm
        else:
            self._set_g_pas(Rin, [self.patch])


class ModelWholeCellPatchCaSingleChan(ModelWholeCellPatchCa):
    def __init__(self, path_to_mods: str, channel_name: str, ion_name: str,
                 external_conc, internal_conc=None,
                 E_rev={}, gbar_name="gbar", temp=22, recompile=True, L=10,
                 diam=10, Ra=100,
                 ljp=0, cvode=True,  Rin=200e6, Rm=20000, cap=None, cm=1,
                 v_rest=-65, gbar_value=0.001,
                 t_decay=20,
                 buffer_capacity=18,
                 membrane_shell_width=memb_shell_width):
        channel_names = [channel_name]
        ion_names = [ion_name]
        
        if external_conc is not None:
            exc = {ion_name: external_conc}
        else:
            raise SystemError("Ca external conc needs to be specified")
        inc = {}
        if internal_conc is not None:
            inc[ion_name] = internal_conc
        gbarn = {channel_name: gbar_name}
        gbarv = {channel_name: gbar_value}
        super(ModelWholeCellPatchCaSingleChan, self).__init__(path_to_mods,
                                                              channel_names,
                                                              ion_names,
                                                              external_conc=exc,
                                                              internal_conc=inc,
                                                              E_rev=E_rev,
                                                              gbar_names=gbarn,
                                                              gbar_values=gbarv,
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


class ModelOocyte(ModelWholeCellPatch):
    #parameters from  PMID: 20737886 DOI: 10.1016/0012-1606(81)90417-6 
    def __init__(self, path_to_mods, channel_names: list, ion_names: list,
                 external_conc={}, internal_conc={}, gbar_names={},
                 temp=22, recompile=True, ljp=0, cvode=True,
                 E_rev={}, gbar_values={}):
        super(ModelOocyte, self).__init__(path_to_mods, channel_names,
                                          ion_names,
                                          external_conc=external_conc,
                                          internal_conc=internal_conc,
                                          gbar_names=gbar_names,
                                          temp=22, recompile=True, L=1.3e3,
                                          diam=1.3e3, Ra=100,
                                          ljp=0, cvode=True,
                                          Rin=1.86e6, Rm=10000,
                                          cap=None, cm=12,
                                          v_rest=-50, E_rev=E_rev,
                                          gbar_values=gbar_values)


class ModelOocyteCa(ModelWholeCellPatchCa):
    #parameters from  PMID: 20737886 DOI: 10.1016/0012-1606(81)90417-6
    # (electric)
    #Ca params from https://doi.org/10.1016/j.ydbio.2005.10.034
    def __init__(self, path_to_mods, channel_names: list, ion_names: list,
                 external_conc={}, internal_conc={}, E_rev={}, gbar_names={},
                 temp=22, recompile=True, ljp=0, cvode=True,
                 gbar_values={}):
        super(ModelOocyteCa, self).__init__(path_to_mods, channel_names,
                                            ion_names,
                                            external_conc=external_conc,
                                            internal_conc=internal_conc,
                                            E_rev=E_rev, gbar_names=gbar_names,
                                            temp=22, recompile=True,
                                            L=1.3e3, diam=1.3e3, Ra=100,
                                            ljp=0, cvode=True,
                                            Rin=1.86e6, Rm=20000,
                                            cap=None, cm=12,
                                            v_rest=-50,
                                            gbar_values=gbar_values,
                                            t_decay=8e3)


class ModelGiantExcisedPatch(ModelWholeCellPatch):
    def __init__(self, path_to_mods, channel_names: list, ion_names: list,
                 external_conc={}, internal_conc={}, gbar_names={},
                 temp=22, recompile=True, ljp=0, cvode=True,
                 E_rev={}, gbar_values={}):
        super(ModelGiantExcisedPatch, self).__init__(path_to_mods,
                                                     channel_names,
                                                     ion_names,
                                                     external_conc=external_conc,
                                                     internal_conc=internal_conc,
                                                     gbar_names=gbar_names,
                                                     temp=temp, recompile=True,
                                                     L=15,
                                                     diam=15, Ra=100,
                                                     ljp=0, cvode=True,
                                                     Rin=None, Rm=400000,
                                                     cap=12e-3, cm=12,
                                                     v_rest=-50, E_rev=E_rev,
                                                     gbar_values=gbar_values)
        self.vclamp.rs = 10e3


class ModelGiantExcisedPatchCa(ModelWholeCellPatchCa):
    def __init__(self, path_to_mods, channel_names: list, ion_names: list,
                 external_conc={}, internal_conc={}, gbar_names={},
                 temp=22, recompile=True, ljp=0, cvode=True,
                 E_rev={}, gbar_values={}):
        super(ModelGiantExcisedPatchCa, self). __init__(path_to_mods,
                                                        channel_names,
                                                        ion_names,
                                                        external_conc=external_conc,
                                                        internal_conc=internal_conc,
                                                        E_rev=E_rev,
                                                        gbar_names=gbar_names,
                                                        temp=temp,
                                                        recompile=recompile,
                                                        ljp=ljp,
                                                        cvode=cvode,
                                                        L=15, diam=15, Ra=100,
                                                        Rin=5e9, cap=None,
                                                        Rm=400000, cm=12e-3,
                                                        v_rest=-50,
                                                        gbar_values=gbar_values,
                                                        t_decay=8, 
                                                        buffer_capacity=10,
                                                        membrane_shell_width=0.1)

        self.vclamp.rs = 10e3


class ModelCellAttachedPatch(ModelWholeCellPatch):
    def __init__(self, path_to_mods, channel_names: list, ion_names: list,
                 external_conc={}, internal_conc={}, gbar_names={},
                 temp=22, recompile=True, ljp=0, cvode=True,
                 E_rev={}, gbar_values={}, Rin=5e9, v_rest=-65):
        super(ModelWholeCellPatch, self).__init__(path_to_mods,
                                                  channel_names,
                                                  ion_names,
                                                  external_conc=external_conc,
                                                  internal_conc=internal_conc,
                                                  gbar_names=gbar_names,
                                                  temp=temp, recompile=True,
                                                  L=3,
                                                  diam=3, Ra=10,
                                                  ljp=0, cvode=True,
                                                  Rin=Rin, Rm=20000,
                                                  cap=None, cm=1,
                                                  v_rest=-65, E_rev=E_rev,
                                                  gbar_values=gbar_values)
        self.vclamp.rs = 10e3


class ModelCellAttachedPatchCa(ModelWholeCellPatchCa):
    def __init__(self, path_to_mods, channel_names: list, ion_names: list,
                 external_conc={}, gbar_names={},
                 temp=22, recompile=True, ljp=0, cvode=True,
                 E_rev={}, gbar_values={}, Rin=5e9, v_rest=-65, t_decay=100,
                 buffer_capacity=20):
        super(ModelCellAttachedPatchCa, self). __init__(path_to_mods,
                                                        channel_names,
                                                        ion_names,
                                                        external_conc=external_conc,
                                                        internal_conc=internal_conc,
                                                        E_rev=E_rev,
                                                        gbar_names=gbar_names,
                                                        temp=temp,
                                                        recompile=recompile,
                                                        ljp=ljp,
                                                        cvode=cvode,
                                                        L=15, diam=15, Ra=100,
                                                        Rin=Rin, cap=None,
                                                        Rm=400, cm=1,
                                                        v_rest=v_rest,
                                                        gbar_values=gbar_values,
                                                        t_decay=t_decay, 
                                                        buffer_capacity=buffer_capacity,
                                                        membrane_shell_width=0.1)

        self.vclamp.rs = 10e3
