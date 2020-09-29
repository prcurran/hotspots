from ccdc import io
import os
from hotspots.pharmacophore_extension import LigandPharmacophoreModel



class PharmacophoreExtensionTD:
    @staticmethod
    def create_multiple_pharmacophore_files():

        # setting output directory
        out_dir = "testdata/pharmacophore_extension/PharmacophoreModel"
        if not os.path.exists(out_dir):
            os.mkdir(out_dir)

        # create files
        mols = io.CrystalReader("testdata/pharmacophore_extension/provided_data/test_overlay.mol2")
        for mol in mols:
            lp = LigandPharmacophoreModel()
            lp.feature_definitions = ["ring_planar_projected"]
            lp.detect_from_ligand(mol)
            for feature in lp.detected_features:
                lp.add_feature(feature)

            lp.write(os.path.join(out_dir, f"{mol.identifier}.cm"))


if __name__ == "__main__":
    ctd = PharmacophoreExtensionTD()
    ctd.create_multiple_pharmacophore_files()