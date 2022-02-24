# SPDX-License-Identifier: LGPL-2.1-or-later
from .patch_models import ModelWholeCellPatch
from .patch_models import ModelWholeCellPatchCa
from .patch_models import ModelWholeCellPatchSingleChan
from .patch_models import ModelWholeCellPatchCaSingleChan
from .patch_models import ModelGiantExcisedPatch
from .patch_models import ModelGiantExcisedPatchCa
from .patch_models import ModelOocytePatch
from .patch_models import ModelOocytePatchCa


import os
loc = os.path.dirname(os.path.abspath(__file__))
data_path = os.path.join(loc, 'demo_CA1')
mechanisms_path = os.path.join(loc, 'mechanisms')
