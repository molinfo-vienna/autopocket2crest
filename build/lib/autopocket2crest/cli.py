import argparse
from .pocket_preparation import run_pipeline

def main():
    parser = argparse.ArgumentParser(
        description="AutoPocket2CREST — Automated protein–ligand pocket preparation for CREST."
    )
    parser.add_argument("protein_file", help="Protein PDB file")
    parser.add_argument("ligand_file", help="Ligand MOL2 file")
    parser.add_argument("pdbid", help="PDB ID (used for directory and naming)")
    parser.add_argument("--no-crest", action="store_true", help="Skip running CREST (prep only)")
    args = parser.parse_args()

    run_pipeline(args.protein_file, args.ligand_file, args.pdbid, run_crest=not args.no_crest)
