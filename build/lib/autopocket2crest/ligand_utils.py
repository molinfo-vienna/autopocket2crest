def get_ligand_name(mol2_file, pdb_file):
    """Detect ligand residue name using MOL2 file."""

    ligand_name = None
    with open(mol2_file, 'r') as f:
        lines = f.readlines()
        next_line = False
        for line in lines:
            if next_line:
                ligand_name = line.split()[1]
                break
            elif line.startswith("@<TRIPOS>SUBSTRUCTURE"):
                next_line = True
    
    if not ligand_name or ligand_name == "UNNAMED":
        with open(mol2_file, 'r') as mol:
            for line in mol:
                if next_line:
                    parts = line.split()
                    if len(parts) > 7:
                        ligand_name = parts[7]
                        break
                elif line.startswith("@<TRIPOS>ATOM"):
                    next_line = True
 
    
    if not ligand_name:
        print(f"Warning: No ligand name found for {mol2_file}, using 'UNNAMED'")
        ligand_name = "UNNAMED"

    return ligand_name
