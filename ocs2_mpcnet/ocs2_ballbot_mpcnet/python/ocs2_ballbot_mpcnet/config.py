###############################################################################
# Copyright (c) 2022, Farbod Farshidian. All rights reserved.
#
# Redistribution and use in source and binary forms, with or without
# modification, are permitted provided that the following conditions are met:
#
#  * Redistributions of source code must retain the above copyright notice, this
#   list of conditions and the following disclaimer.
#
#  * Redistributions in binary form must reproduce the above copyright notice,
#   this list of conditions and the following disclaimer in the documentation
#   and/or other materials provided with the distribution.
#
#  * Neither the name of the copyright holder nor the names of its
#   contributors may be used to endorse or promote products derived from
#   this software without specific prior written permission.
#
# THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
# AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
# IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
# DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE
# FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL
# DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR
# SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER
# CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY,
# OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
# OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
###############################################################################

"""Ballbot configuration variables.

Sets robot-specific configuration variables for ballbot.
"""

from ocs2_mpcnet_core import config

#
# config
#

# data type for tensor elements
DTYPE = config.DTYPE

# device on which tensors will be allocated
DEVICE = config.DEVICE

#
# ballbot_config
#

# name of the robot
NAME = "ballbot"

# (generalized) time dimension
TIME_DIM = 1

# state dimension
STATE_DIM = 10

# input dimension
INPUT_DIM = 3

# target trajectories state dimension
TARGET_STATE_DIM = STATE_DIM

# target trajectories input dimension
TARGET_INPUT_DIM = INPUT_DIM

# observation dimension
OBSERVATION_DIM = STATE_DIM

# action dimension
ACTION_DIM = INPUT_DIM

# observation scaling
# fmt: off
OBSERVATION_SCALING = [
    1.0, 1.0, 1.0, 1.0, 1.0,  # pose
    1.0, 1.0, 1.0, 1.0, 1.0   # twist
]
# fmt: on

# action scaling
ACTION_SCALING = [1.0, 1.0, 1.0]

# input cost for behavioral cloning
R = [2.0, 2.0, 2.0]
