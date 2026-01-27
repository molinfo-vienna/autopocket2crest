def filter_pdb_by_altloc(input_file, output_file, keep_altloc='A'):
    """Keep only the desired alternate location (altLoc) in a PDB file."""
    with open(input_file, 'r') as fin, open(output_file, 'w') as fout:
        for line in fin:
            if line.startswith(('ATOM', 'HETATM')):
                altloc = line[16]
                if altloc == ' ' or altloc == keep_altloc:
                    fout.write(line)
            else:
                fout.write(line)

def cut_pocket(ligand, ligand_resname, u):

    if len(ligand) == 0:
        print(f"Warning: No ligand found with residue name {ligand_resname}. Exiting.")
        exit(1)
    elif len(ligand) > 120:
        print(f"Warning: More than 120 atoms found in the ligand. This might be an error. Exiting.")
        exit(1)
    else:
        print(f"Ligand found with {len(ligand)} atoms.")

    print("Selecting pocket atoms...")

    pocket = u.select_atoms("around 3 group ligand_group", ligand_group=ligand)
    print(f"Initial pocket atoms (first shell): {len(pocket)}")

    i = 0
    while len(pocket) < 70:
        print("No atoms found in the pocket. Increasing the selection radius...")
        # Increase the radius by 0.5 Å and try again
        i += 0.5
        pocket = u.select_atoms(f"around {3 + i} group ligand_group", ligand_group=ligand)
        if i > 50:  # Prevent infinite loop
            print("Warning: No atoms found in the pocket after 50 attempts. Exiting.")
            exit(1)
        else:
            print(f"Found pocket atoms after {i*2} iterations with radius {3 + i} Å")
            print(f"New Pocket atoms (first shell): {len(pocket)}")

    pocket = pocket.residues.atoms
    # === Expand pocket 2.6 Å, exclude ligand again ===
    pocket_extended = u.select_atoms(
        "group pocket_group or (around 2.6 group pocket_group and not group ligand_group)",
        pocket_group=pocket,
        ligand_group=ligand
    )

    print(f"Extended pocket atoms: {len(pocket_extended)}")


    # === Write pocket structure ===
    pocket_extended = pocket_extended.select_atoms(f"not resname {ligand_resname} and not resname HOH")

    isolated_atoms = []
    for i, atom in enumerate(pocket_extended.atoms):
        atom_group = pocket_extended.atoms[[i]]  # Use local index within AtomGroup
        test_group = pocket_extended.select_atoms(f"around 1.9 group atom", atom=atom_group)
        if len(test_group) == 0:
            print(f"Warning: Atom {atom} is isolated in the pocket. Removing it.")
            isolated_atoms.append(atom)
    for atom2 in isolated_atoms:
            pocket_extended = pocket_extended.select_atoms(f"not id {atom2.id}")

    pocket_extended.write("test_pocket_extended.pdb")
    print(f"Pocket size: {len(pocket_extended)} atoms.")
    
def fix_pdb_elements(input_pdb, output_pdb):
    """
    Fix missing element identifiers in a PDB file by inferring them
    from the atom name column. Lines that already have an element
    identifier are left unchanged.

    Parameters
    ----------
    input_pdb : str
        Path to the input PDB file
    output_pdb : str
        Path to the output PDB file with fixed element columns
    """

    def infer_element(atom_name):
        # Remove digits and spaces
        name = ''.join(c for c in atom_name if c.isalpha())

        if not name:
            return "  "

        # Two-letter elements have a lowercase second letter in atom name
        if len(name) >= 2 and name[1].islower():
            return name[:2].capitalize()
        else:
            return name[0].upper()

    with open(input_pdb, "r") as fin, open(output_pdb, "w") as fout:
        for line in fin:
            if line.startswith(("ATOM", "HETATM")):
                # Make sure the line is long enough to contain element field
                padded = line.rstrip("\n").ljust(78)

                # Element field is columns 77–78 (0-based 76:78)
                existing_element = padded[76:78].strip()

                if existing_element:
                    # Element already present → write original line unchanged
                    fout.write(line)
                else:
                    atom_name = padded[12:16]
                    element = infer_element(atom_name)

                    fixed_line = padded[:76] + element.rjust(2) + padded[78:]
                    fout.write(fixed_line + "\n")
            else:
                fout.write(line)
