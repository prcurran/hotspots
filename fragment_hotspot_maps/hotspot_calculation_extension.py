'''
The main class of the :mod:`fragment_hotspot_maps.hotspot_calculation_extension.Hotspots`.

This is an internal extension of :class:`ccdc.grid.Grid` that adds in potential new features
for review and development. These features are based on non-public dependancies
'''
#############################################################################

from ccdc_rp import cavity

import hotspot_calculation

class Hotspots(hotspot_calculation.Hotspots()):

    def from_cavity(self):
        hr = 0
        return hr

    def _get_cavity_centroid(self, cav):
        """
        Returns the centroid of a cavity object

        :param cav:
        :return:
        """

        x_coords = []
        y_coords = []
        z_coords = []

        for feat in cav.features:
            feature_coords = feat.coordinates
            x_coords.append(feature_coords[0])
            y_coords.append(feature_coords[1])
            z_coords.append(feature_coords[2])

        x_avg = round(np.mean(x_coords))
        y_avg = round(np.mean(y_coords))
        z_avg = round(np.mean(z_coords))

        return x_avg, y_avg, z_avg

hotspot_calculation.Hotspots = Hotspots

# def _score_cavity_features(self, cav):
#     '''
#     Assign a score to each cavity feature, based on the feature's ideal interaction point
#
#     :param cav: a :class:`ccdc.Cavity.cavity` instance
#     :return: a dictionary {score:[list of features], where each item in list of features is a tuple of (atom, feature, atom_type)
#     '''
#
#     '''Assigns a score to each cavity feature. Donor scores are assigned to polar hydrogens, rather than the heavy
#         atom. Returns a dictionary of {score:[features]}'''
#
#     for feat in cav.features:
#         atom_type_dict = {'pi': 'apolar', 'aromatic': 'apolar', 'aliphatic': 'apolar', 'donor': 'acceptor',
#                           'acceptor': 'donor'}
#         feat_coords = feat.coordinates
#         atom = self._get_feat_atom(feat_coords)
#         if atom is None:
#             continue
#         ideal_vector = feat.protein_vector
#         coords = self._get_ideal_coordinates(ideal_vector, feat_coords, 3)
#         if feat.type == 'donor_acceptor':
#             d_score = self._get_near_score(coords, 'acceptor', tolerance=4)
#             a_score = self._get_near_score(coords, 'donor', tolerance=4)
#             self._update_score_dic(self.features_by_score, atom, feat, a_score, 'acceptor')
#             for n in atom.neighbours:
#                 if n.atomic_number == 1:
#                     self._update_score_dic(self.features_by_score, n, feat, d_score, 'donor')
#
#         elif atom_type_dict[feat.type] == 'acceptor':
#             score = self._get_near_score(coords, atom_type_dict[feat.type], tolerance=3)
#             for n in atom.neighbours:
#                 if n.atomic_number == 1:
#                     self._update_score_dic(self.features_by_score, n, feat, score, 'donor')
#         elif atom_type_dict[feat.type] == 'donor':
#             score = self._get_near_score(coords, atom_type_dict[feat.type], tolerance=3)
#             self._update_score_dic(self.features_by_score, atom, feat, score, 'acceptor')
#
#return self.features_by_score