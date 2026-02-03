import os, time, MDAnalysis as mda
from rdkit import Chem
from .pdb_utils import filter_pdb_by_altloc, cut_pocket, fix_pdb_elements
from .ligand_utils import get_ligand_name
from .remove_unbonded import remove_unbonded
from .crest_interface import generate_constraints, run_crest
from .transfer_pdb_info import transfer_pdb_info
from .cleanup import cleanup_temp_files

def run_pipeline(protein_file, ligand_file, outdir, run_crest_bool=True, base_dir=".", temp="310", lvl_of_theory="gfnff", extra_crest_args="-squick"):
    """Full AutoPocket2CREST pipeline."""
    print("Starting AutoPocket2CREST pipeline...")    
    start = time.time()
    os.makedirs(os.path.join(base_dir, outdir), exist_ok=True)
    os.chdir(os.path.join(base_dir, outdir))

    cwd = os.getcwd()
    print(f"Current working directory: {cwd}")
    
    fix_pdb_elements(f"{base_dir}/{protein_file}", f"{cwd}/pre_prepared.pdb")
    filter_pdb_by_altloc(f"{cwd}/pre_prepared.pdb", "prepared.pdb")
    ligand_resname = get_ligand_name(f"{base_dir}/{ligand_file}", f"{base_dir}/{protein_file}")
    print(f"Ligand identified: {ligand_resname}")

    u_prot = mda.Universe("prepared.pdb")
    u_lig = mda.Universe(f"{base_dir}/{ligand_file}")
    u = mda.Merge(u_prot.select_atoms("protein"), u_lig.atoms)
    ligand = u.select_atoms(f"resname {ligand_resname}")

    print("Cutting pocket around the ligand...")
    cut_pocket(ligand, ligand_resname, u)
    os.system(f"obabel -ipdb {cwd}/test_pocket_extended.pdb -opdb -O {cwd}/test_pocket_extended.pdb -d")
    remove_unbonded(f"{cwd}/test_pocket_extended.pdb", f"{cwd}/test_pocket_extended_clean.pdb")

    # Hydrogenation (Open Babel)
    print("Adding hydrogens to pocket...")
    os.system(f"obabel -ipdb {cwd}/test_pocket_extended_clean.pdb -opdb -O {cwd}/test_pocket_extended_h.pdb -p 7.4")

    # Merge ligand + pocket_h
    print("Merging ligand and hydrogenated pocket...")
    u_pocket = mda.Universe(f"{cwd}/test_pocket_extended_h.pdb")
    full = mda.Merge(ligand, u_pocket.atoms)
    full.atoms.write(f"{cwd}/test_pocket_extended_h_fixed.pdb")

    # Calculate charge
    print("Calculating formal charge of the system...")
    mol = Chem.MolFromPDBFile(f"{cwd}/test_pocket_extended_h_fixed.pdb", sanitize=False, removeHs=False)
    charge = Chem.GetFormalCharge(mol)
    print("Formal charge:", charge)

    if run_crest_bool:
        print("Running CREST conformer search...")
        # Select backbone atoms from this definitive final structure
        u_final = mda.Universe(f"{cwd}/test_pocket_extended_h_fixed.pdb")
        backbone_sel = u_final.select_atoms(f"not resname {ligand_resname}")
        backbone = (backbone_sel.indices + 1).tolist()
        constraint_file = generate_constraints(f"{cwd}/test_pocket_extended_h_fixed.pdb", backbone)
        run_crest(f"{cwd}/test_pocket_extended_h_fixed.pdb", constraint_file, charge, temp=temp, lvl_of_theory=lvl_of_theory, extra_crest_args=extra_crest_args)
        os.system("obabel -ixyz crest_conformers.xyz -opdb -O crest_conformers.pdb")

        transfer_pdb_info(f"{cwd}/test_pocket_extended_h_fixed.pdb", f"{cwd}/crest_conformers.pdb", f"{cwd}/crest_conformers_updated.pdb")
        print("CREST conformer search complete.")

    cleanup_temp_files([
    f"{cwd}/test_pocket_extended.pdb",
    f"{cwd}/test_pocket_extended_h.pdb",
    f"{cwd}/pre_prepared.pdb",
    f"{cwd}/prepared.pdb",
    f"{cwd}/.CHRG",
    f"{cwd}/.xcontrol.sample"
    ])
    print("AutoPocket2CREST complete in", round((time.time()-start)/3600, 2), "hours.")
