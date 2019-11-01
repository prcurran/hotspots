"""
Tutorial script for using the ensemble and selectivity maps
"""
from hotspots.hs_ensembles import EnsembleResult, SelectivityResult
from hotspots.hs_io import HotspotReader, HotspotWriter
from pathlib import Path

# Use a list of precalculated hotspot results to make ensmeble maps for p38alpha and ERK2 kinases.

# Set the on and off targets and load the pre-calculated hotspots results
on_target = "p38alpha"
off_target = "ERK2"

ensemble_settings = EnsembleResult.Settings()
ensemble_settings.combine_mode = 'median'
ensemble_settings.apolar_frequency_threshold = 0.0
ensemble_settings.polar_frequency_threshold = 20.0

ensemble_results = {}

for t in [on_target, off_target]:
    target_paths = list(Path(f"{t}").glob('*/fullsize_hotspots_3000/out.zip'))

    # Compile into lists of hotspot results
    ensemble = EnsembleResult(hs_results_list=[HotspotReader(p).read() for p in target_paths],
                                    ensemble_id=on_target,
                                    settings=ensemble_settings)

    ensemble.make_ensemble_maps(save_grid_ensembles=False)
    ensemble_hs_result = ensemble.ensemble_hotspot_result
    with HotspotWriter(f"ensemble_{t}", grid_extension=".ccp4", zip_results=False) as w:
        w.write(ensemble_hs_result)
    ensemble_results[t] = ensemble_hs_result


selectivity_settings = SelectivityResult.Settings()
selectivity_settings.cluster_distance_cutoff = 3.0
selectivity_settings.minimal_cluster_score = 10.0

p38alpha_over_ERK2 = SelectivityResult(target_result=ensemble_results[on_target],
                                       other_result=ensemble_results[off_target])
p38alpha_over_ERK2.make_selectivity_maps()

with HotspotWriter(f"selectivity_{on_target}_over_{off_target}",grid_extension=".ccp4", zip_results=False) as w:
    w.write(p38alpha_over_ERK2.selectivity_result)