# SPDX-License-Identifier: LGPL-2.1-or-later
from .model_patch import ModelPatchNernst
from .model_patch import ModelWholeCellPatchNernst
from .model_patch import ModelPatchConcentration
from .model_patch import ModelWholeCellPatchConcentration
import os
loc = os.path.dirname(os.path.abspath(__file__))
data_path = os.path.join(loc, 'demo_CA1')
mechanisms_path = os.path.join(loc, 'mechanisms')
