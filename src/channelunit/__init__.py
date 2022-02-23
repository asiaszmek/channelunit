# SPDX-License-Identifier: LGPL-2.1-or-later
from .model_patches import ModelWholeCellPatch
from .model_patches import ModelWholeCellPatchCa
from .model_patches import ModelWholeCellPatchSingleChan
from .model_patches import ModelWholeCellPatchCaSingleChan
from .model_patches import ModelGiantExcisedPatch
from .model_patches import ModelGiantExcisedPatchCa
from .model_patches import ModelOocytePatch
from .model_patches import ModelOocytePatchCa


import os
loc = os.path.dirname(os.path.abspath(__file__))
data_path = os.path.join(loc, 'demo_CA1')
mechanisms_path = os.path.join(loc, 'mechanisms')
