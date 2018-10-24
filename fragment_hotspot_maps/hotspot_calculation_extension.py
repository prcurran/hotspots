'''
The main class of the :mod:`fragment_hotspot_maps.hotspot_calculation_extension.Hotspots`.

This is an internal extension of :class:`ccdc.grid.Grid` that adds in potential new features
for review and development. These features are based on non-public dependancies
'''
#############################################################################

import hotspot_calculation


# def _get_cavity_centroid(self, cav):
#     """
#     Returns the centroid of a cavity object
#
#     :param cav:
#     :return:
#     """
#
#     x_coords = []
#     y_coords = []
#     z_coords = []
#
#     for feat in cav.features:
#         feature_coords = feat.coordinates
#         x_coords.append(feature_coords[0])
#         y_coords.append(feature_coords[1])
#         z_coords.append(feature_coords[2])
#
#     x_avg = round(np.mean(x_coords))
#     y_avg = round(np.mean(y_coords))
#     z_avg = round(np.mean(z_coords))
#
#     return x_avg, y_avg, z_avg

class Hotspots(hotspot_calculation.Hotspots()):

    def from_cavity(self):
        hr = 0
        return hr

hotspot_calculation.Hotspots = Hotspots
