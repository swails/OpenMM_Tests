"""
Houses all of the "extra" test cases to run for OpenMM
"""

__all__ = ['testhewl', 'available_platforms']

import testhewl

# Construct the list of available platforms
available_platforms = []
from simtk.openmm import Platform
for i in range(Platform.getNumPlatforms()):
    available_platforms.append(Platform.getPlatform(i).getName())

# Cleanup the namespace a bit
del i, Platform
