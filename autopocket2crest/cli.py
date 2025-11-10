import argparse
import os
from .pocket_preparation import run_pipeline

def main():
    parser = argparse.ArgumentParser(
        description="AutoPocket2CREST — Automated protein–ligand pocket preparation for CREST."
    )
    parser.add_argument("protein_file", help="Protein PDB file")
    parser.add_argument("ligand_file", help="Ligand MOL2 file")
    parser.add_argument("pdbid", help="PDB ID (used for directory and naming)")
    parser.add_argument("--no-crest", action="store_true", help="Skip running CREST (prep only)")
    parser.add_argument("temp", nargs='?', default="310", help="Temperature for CREST (default: 310 K)")
    parser.add_argument("lvl_of_theory",nargs='?', default="gfnff", help="Level of theory for CREST (default: gfnff)")
    parser.add_argument("extra_crest_args", nargs='?', default="-squick", help="Additional arguments to pass to CREST (default: -squick)")

    args = parser.parse_args()
    directory = os.getcwd()
    
    run_pipeline(args.protein_file, args.ligand_file, args.pdbid, run_crest_bool=not args.no_crest, base_dir=directory,temp=args.temp, lvl_of_theory=args.lvl_of_theory, extra_crest_args=args.extra_crest_args)
    print(f"AutoPocket2CREST pipeline completed. Results are in the directory: {os.path.join(directory, args.pdbid)}")
