# SPDX-License-Identifier: LGPL-2.1-or-later
from .model_patch import ModelWholeCellPatch
from .model_patch import ModelWholeCellPatchCa
from .model_patch import ModelWholeCellPatchSingleChan
from .model_patch import ModelWholeCellPatchCaSingleChan
import os
loc = os.path.dirname(os.path.abspath(__file__))
data_path = os.path.join(loc, 'demo_CA1')
mechanisms_path = os.path.join(loc, 'mechanisms')
