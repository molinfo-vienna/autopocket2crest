import os, time, subprocess, MDAnalysis as mda
from rdkit import Chem
from .pdb_utils import filter_pdb_by_altloc, cut_pocket, fix_pdb_elements
from .ligand_utils import get_ligand_name
from .remove_unbonded import remove_unbonded
from .crest_interface import generate_constraints, run_crest
from .transfer_pdb_info import transfer_pdb_info
from .cleanup import cleanup_temp_files

def run_pipeline(
    protein_file,
    ligand_file,
    outdir,
    run_crest_bool=True,
    base_dir=".",
    temp="310",
    lvl_of_theory="gfnff",
    extra_crest_args="-squick",
):
    """Full AutoPocket2CREST pipeline."""

    print("Starting AutoPocket2CREST pipeline...")
    start = time.time()

    os.makedirs(
        os.path.join(base_dir, outdir),
        exist_ok=True,
    )

    os.chdir(
        os.path.join(base_dir, outdir)
    )

    cwd = os.getcwd()

    print(f"Current working directory: {cwd}")

    # ------------------------------------------------------------
    # Prepare input structures
    # ------------------------------------------------------------

    fix_pdb_elements(
        f"{base_dir}/{protein_file}",
        f"{cwd}/pre_prepared.pdb",
    )

    filter_pdb_by_altloc(
        f"{cwd}/pre_prepared.pdb",
        "prepared.pdb",
    )

    ligand_resname = get_ligand_name(
        f"{base_dir}/{ligand_file}",
        f"{base_dir}/{protein_file}",
    )

    print(f"Ligand identified: {ligand_resname}")

    u_prot = mda.Universe("prepared.pdb")
    u_lig = mda.Universe(
        f"{base_dir}/{ligand_file}"
    )

    u = mda.Merge(
        u_prot.select_atoms("protein"),
        u_lig.atoms,
    )

    ligand = u.select_atoms(
        f"resname {ligand_resname}"
    )

    # ------------------------------------------------------------
    # Cut pocket
    # ------------------------------------------------------------

    print("Cutting pocket around the ligand...")

    cut_pocket(
        ligand,
        ligand_resname,
        u,
    )

    result = subprocess.run(
        [
            "obabel",
            "-ipdb",
            f"{cwd}/test_pocket_extended.pdb",
            "-opdb",
            "-O",
            f"{cwd}/test_pocket_extended.pdb",
            "-d",
        ],
        capture_output=True,
        text=True,
    )

    if result.returncode != 0:
        raise RuntimeError(
            "Open Babel failed while preparing the pocket:\n"
            f"{result.stderr}"
        )

    remove_unbonded(
        f"{cwd}/test_pocket_extended.pdb",
        f"{cwd}/test_pocket_extended_clean.pdb",
    )

    # ------------------------------------------------------------
    # Hydrogenation
    # ------------------------------------------------------------

    print("Adding hydrogens to pocket...")

    result = subprocess.run(
        [
            "obabel",
            "-ipdb",
            f"{cwd}/test_pocket_extended_clean.pdb",
            "-opdb",
            "-O",
            f"{cwd}/test_pocket_extended_h.pdb",
            "-p",
            "7.4",
        ],
        capture_output=True,
        text=True,
    )

    if result.returncode != 0:
        raise RuntimeError(
            "Open Babel failed during hydrogenation:\n"
            f"{result.stderr}"
        )

    # ------------------------------------------------------------
    # Merge ligand + hydrogenated pocket
    # ------------------------------------------------------------

    print("Merging ligand and hydrogenated pocket...")

    u_pocket = mda.Universe(
        f"{cwd}/test_pocket_extended_h.pdb"
    )

    full = mda.Merge(
        ligand,
        u_pocket.atoms,
    )

    full.atoms.write(
        f"{cwd}/test_pocket_extended_h_fixed.pdb"
    )

    # ------------------------------------------------------------
    # Calculate charge
    # ------------------------------------------------------------

    print("Calculating formal charge of the system...")

    mol = Chem.MolFromPDBFile(
        f"{cwd}/test_pocket_extended_h_fixed.pdb",
        sanitize=False,
        removeHs=False,
    )

    if mol is None:
        raise RuntimeError(
            "RDKit could not read test_pocket_extended_h_fixed.pdb"
        )

    charge = Chem.GetFormalCharge(mol)

    print("Formal charge:", charge)

    print(
        "Preparation complete in",
        round(time.time() - start, 2),
        "seconds",
    )

    # ============================================================
    # CREST
    # ============================================================

    if run_crest_bool:

        print("Running CREST conformer search...")

        # --------------------------------------------------------
        # Remove stale CREST output
        # --------------------------------------------------------

        for filename in [
            f"{cwd}/crest_conformers.xyz",
            f"{cwd}/crest_conformers.pdb",
            f"{cwd}/crest_conformers_updated.pdb",
            f"{cwd}/crest.out",
        ]:
            if os.path.exists(filename):
                os.remove(filename)
                print(
                    f"Removed stale CREST output: {filename}"
                )

        # --------------------------------------------------------
        # Select constrained backbone atoms
        # --------------------------------------------------------

        u_final = mda.Universe(
            f"{cwd}/test_pocket_extended_h_fixed.pdb"
        )

        backbone_sel = u_final.select_atoms(
            f"not resname {ligand_resname}"
        )

        backbone = (
            backbone_sel.indices + 1
        ).tolist()

        # --------------------------------------------------------
        # Generate CREST constraint file
        # --------------------------------------------------------

        constraint_file = generate_constraints(
            f"{cwd}/test_pocket_extended_h_fixed.pdb",
            backbone,
        )

        # --------------------------------------------------------
        # Run CREST
        # --------------------------------------------------------

        run_crest(
            f"{cwd}/test_pocket_extended_h_fixed.pdb",
            constraint_file,
            charge,
            temp=temp,
            lvl_of_theory=lvl_of_theory,
            extra_crest_args=extra_crest_args,
        )

        # --------------------------------------------------------
        # Verify CREST output
        # --------------------------------------------------------

        crest_xyz = (
            f"{cwd}/crest_conformers.xyz"
        )

        if not os.path.isfile(crest_xyz):
            raise RuntimeError(
                "CREST completed without producing "
                "crest_conformers.xyz."
            )

        # --------------------------------------------------------
        # Convert XYZ -> PDB
        # --------------------------------------------------------

        print(
            "Converting CREST conformers to PDB..."
        )

        crest_pdb = (
            f"{cwd}/crest_conformers.pdb"
        )

        result = subprocess.run(
            [
                "obabel",
                "-ixyz",
                crest_xyz,
                "-opdb",
                "-O",
                crest_pdb,
            ],
            capture_output=True,
            text=True,
        )

        if result.returncode != 0:
            raise RuntimeError(
                "Open Babel failed while converting "
                "CREST output to PDB:\n"
                f"{result.stderr}"
            )

        # --------------------------------------------------------
        # Transfer PDB information
        # --------------------------------------------------------

        print(
            "Transferring PDB atom/residue information..."
        )

        transfer_pdb_info(
            f"{cwd}/test_pocket_extended_h_fixed.pdb",
            crest_pdb,
            f"{cwd}/crest_conformers_updated.pdb",
        )

        print(
            "CREST conformer search complete."
        )

    # ============================================================
    # Cleanup
    # ============================================================

    cleanup_temp_files([
        f"{cwd}/test_pocket_extended.pdb",
        f"{cwd}/test_pocket_extended_h.pdb",
        f"{cwd}/pre_prepared.pdb",
        f"{cwd}/prepared.pdb",
        f"{cwd}/.CHRG",
        f"{cwd}/.xcontrol.sample",
    ])

    print(
        "Full AutoPocket2CREST Pipeline complete in",
        round(
            (time.time() - start) / 3600,
            2,
        ),
        "hours.",
    )