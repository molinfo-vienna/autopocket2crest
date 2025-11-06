def remove_unbonded(pdb, output):
    """Remove residues that are not covalently bonded to others."""
    with open(pdb, 'r') as f:
        lines = f.readlines()

    atom_info = {}
    residue_atoms = {}
    conect_map = {}

    # --- Parse atoms and CONECT records ---
    for line in lines:
        if line.startswith(("ATOM", "HETATM")):
            res_name = line[17:20].strip()
            chain_id = line[21].strip()
            res_seq = int(line[22:26].strip())
            res_id = (chain_id, res_name, res_seq)
            atom_serial = int(line[6:11])
            atom_info[atom_serial] = (res_id, line)
            residue_atoms.setdefault(res_id, set()).add(atom_serial)
        elif line.startswith("CONECT"):
            a1 = int(line[6:11])
            bonded = [line[i:i+5].strip() for i in (11, 16, 21, 26)]
            for a2_str in bonded:
                if a2_str:
                    a2 = int(a2_str)
                    conect_map.setdefault(a1, set()).add(a2)
                    conect_map.setdefault(a2, set()).add(a1)

    # --- Identify residues connected to others ---
    residues_to_keep = set()
    for atom_serial, (res_id, _) in atom_info.items():
        for b in conect_map.get(atom_serial, []):
            bonded_res = atom_info.get(b, (None,))[0]
            if bonded_res and bonded_res != res_id:
                residues_to_keep.add(res_id)
                break

    # --- Write output ---
    with open(output, 'w') as f:
        for line in lines:
            if line.startswith(("ATOM", "HETATM")):
                atom_serial = int(line[6:11])
                res_id = atom_info.get(atom_serial, (None,))[0]
                if res_id in residues_to_keep:
                    f.write(line)
            elif line.startswith("CONECT"):
                try:
                    a1 = int(line[6:11])
                    bonded = [int(line[i:i+5]) for i in (11, 16, 21, 26) if line[i:i+5].strip()]
                    all_atoms = [a1] + bonded
                    if all(atom_info.get(a, (None,))[0] in residues_to_keep for a in all_atoms):
                        f.write(line)
                except ValueError:
                    continue
            else:
                f.write(line)

    print(f"Removed unbonded residues. Output: {output}")
    return output
