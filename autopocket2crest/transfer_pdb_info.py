def get_template_fields(template_file):
    data = []
    with open(template_file, 'r') as f:
        for line in f:
            if line.startswith(('ATOM', 'HETATM')):
                data.append((line[12:16], line[17:20], line[22:26]))
    return data

def parse_conformers(conformers_file):
    conformers, current, inside = [], [], False
    with open(conformers_file, 'r') as f:
        for line in f:
            if line.startswith("MODEL"):
                current, inside = [line], True
            elif line.startswith("ENDMDL"):
                current.append(line)
                conformers.append(current)
                inside = False
            elif inside:
                current.append(line)
    if not conformers:
        print("Warning: No models found, treating file as single conformer.")
        with open(conformers_file, 'r') as f:
            conf = ["MODEL        1\n"]
            for line in f:
                if line.startswith(("ATOM","HETATM")):
                    conf.append(line)
            conf.append("ENDMDL\n")
        conformers.append(conf)
    return conformers

def update_model(model_lines, template_data):
    out = []
    for i, line in enumerate(model_lines):
        if line.startswith(('ATOM', 'HETATM')) and i < len(template_data):
            atom_name, res_name, res_num = template_data[i]
            new_line = (
                line[:12] + atom_name + line[16:17] +
                res_name + line[20:22] + res_num + line[26:]
            )
            out.append(new_line)
        else:
            out.append(line)
    return out

def transfer_pdb_info(template_path, conformers_path, output_path):
    """Align conformer residue/atom info with reference template PDB."""
    template_data = get_template_fields(template_path)
    conformers = parse_conformers(conformers_path)
    updated = []
    for model in conformers:
        updated.extend(update_model(model, template_data))
    with open(output_path, 'w') as f:
        f.writelines(updated)
    print(f"Updated conformers written to {output_path}")
    return output_path
