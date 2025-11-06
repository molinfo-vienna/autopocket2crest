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
