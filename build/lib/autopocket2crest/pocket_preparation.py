import os, time, MDAnalysis as mda
from rdkit import Chem
from .pdb_utils import filter_pdb_by_altloc
from .ligand_utils import get_ligand_name
from .remove_unbonded import remove_unbonded
from .crest_interface import generate_constraints, run_crest
from .transfer_pdb_info import transfer_pdb_info
from .cleanup import cleanup_temp_files

def run_pipeline(protein_file, ligand_file, pdbid, run_crest=True):
    """Full AutoPocket2CREST pipeline."""
    start = time.time()
    os.makedirs(pdbid, exist_ok=True)
    os.chdir(pdbid)

    filter_pdb_by_altloc(protein_file, "prepared.pdb")
    ligand_resname = get_ligand_name(ligand_file, protein_file)
    print(f"Ligand identified: {ligand_resname}")

    u_prot = mda.Universe("prepared.pdb")
    u_lig = mda.Universe(ligand_file)
    u = mda.Merge(u_prot.select_atoms("protein"), u_lig.atoms)
    ligand = u.select_atoms(f"resname {ligand_resname}")

    pocket = u.select_atoms("around 3 group ligand_group", ligand_group=ligand)
    while len(pocket) < 70:
        pocket = u.select_atoms(f"around {3 + len(pocket)/50:.1f} group ligand_group", ligand_group=ligand)

    pocket = pocket.residues.atoms
    pocket.write("test_pocket_extended.pdb")

    remove_unbonded("test_pocket_extended.pdb", "test_pocket_extended_clean.pdb")

    # Hydrogenation (Open Babel)
    os.system("obabel -ipdb test_pocket_extended_clean.pdb -opdb -O test_pocket_extended_h.pdb -p 7.4")

    # Merge ligand + pocket_h
    u_pocket = mda.Universe("test_pocket_extended_h.pdb")
    full = mda.Merge(ligand, u_pocket.atoms)
    full.atoms.write("test_pocket_extended_h_fixed.pdb")

    # Calculate charge
    mol = Chem.MolFromPDBFile("test_pocket_extended_h_fixed.pdb", sanitize=False, removeHs=False)
    charge = Chem.GetFormalCharge(mol)
    print("Formal charge:", charge)

    if run_crest:
        backbone = list(range(1, len(full.atoms)//2))  # Simplified
        constraint_file = generate_constraints("test_pocket_extended_h_fixed.pdb", backbone)
        run_crest("test_pocket_extended_h_fixed.xyz", constraint_file, charge)
        transfer_pdb_info("test_pocket_extended_h_fixed.pdb", "crest_conformers.pdb", "crest_conformers_updated.pdb")

    cleanup_temp_files([
    "test_pocket_extended.pdb",
    "test_pocket_extended_h.pdb",
    "prepared.pdb",
    ".CHRG",
    ".xcontrol.sample"
    ])
    print("AutoPocket2CREST complete in", round((time.time()-start)/3600, 2), "hours.")
