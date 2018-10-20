'''
The main class of the :mod:`fragment_hotspot_maps.hotspot_calculation_extension.Hotspots`.

This is an internal extension of :class:`ccdc.grid.Grid` that adds in potential new features
for review and development. These features are based on non-public dependancies
'''
#############################################################################

import hotspot_calculation


class Hotspots(hotspot_calculation.Hotspots()):

    def from_cavity(self):
        hr = 0
        return hr

hotspot_calculation.Hotspots = Hotspots