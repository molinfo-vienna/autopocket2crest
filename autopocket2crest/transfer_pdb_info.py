from collections import defaultdict


def parse_template(template_file):
    """
    Read atom metadata and connectivity from the reference/template PDB.

    The template's existing residue numbers are deliberately ignored.
    """

    atoms = {}
    bonds = defaultdict(set)

    with open(template_file, "r") as f:
        for line in f:

            if line.startswith(("ATOM  ", "HETATM")):
                serial = int(line[6:11])

                atoms[serial] = {
                    "record": line[:6],
                    "atom_name": line[12:16],
                    "altloc": line[16:17],
                    "res_name": line[17:20],
                    "chain": line[21:22],
                    "res_num": line[22:26],
                    "icode": line[26:27],
                    "element": line[76:78],
                    "charge": line[78:80],
                }

            elif line.startswith("CONECT"):
                serial = int(line[6:11])

                for i in range(11, len(line), 5):
                    field = line[i:i+5].strip()

                    if field:
                        other = int(field)

                        # Store the bond in both directions
                        bonds[serial].add(other)
                        bonds[other].add(serial)

    return atoms, bonds


def find_connected_components(atoms, bonds):
    """
    Find molecular components from the CONECT graph.

    Returns a list of lists of atom serial numbers.
    """

    remaining = set(atoms)
    components = []

    while remaining:

        start = min(remaining)
        stack = [start]
        component = set()

        while stack:

            atom = stack.pop()

            if atom in component:
                continue

            component.add(atom)
            remaining.discard(atom)

            for neighbor in bonds.get(atom, set()):
                if neighbor not in component:
                    stack.append(neighbor)

        components.append(sorted(component))

    # Sort components according to their first atom
    components.sort(key=lambda x: x[0])

    return components


def make_residue_mapping(components):
    """
    Assign new residue numbers based on connected components.

    Component 1 -> residue 1
    Component 2 -> residue 2
    etc.
    """

    residue_map = {}

    for residue_number, component in enumerate(components, start=1):

        for atom_serial in component:
            residue_map[atom_serial] = residue_number

    return residue_map


def get_component_names(template_atoms, components):
    """
    Determine a sensible residue name for each component.

    The first component in your template is FMS.
    The other components are currently UNL.
    """

    names = []

    for component in components:

        component_names = {
            template_atoms[serial]["res_name"].strip()
            for serial in component
        }

        # If there is exactly one meaningful residue name, use it.
        if len(component_names) == 1:
            name = next(iter(component_names))
        else:
            name = "UNL"

        names.append(name)

    return names


def parse_conformer_pdb(conformers_file):
    """
    Read the Open Babel multi-model PDB.

    Returns:
        list of models
        each model is a list of lines
    """

    models = []
    current = None

    with open(conformers_file, "r") as f:

        for line in f:

            if line.startswith("MODEL"):
                current = [line]

            elif line.startswith("ENDMDL"):
                if current is not None:
                    current.append(line)
                    models.append(current)
                    current = None

            elif current is not None:
                current.append(line)

    if not models:
        raise ValueError("No MODEL records found in conformers PDB.")

    return models


def update_model(
    model_lines,
    template_atoms,
    residue_map,
    component_names,
):
    """
    Transfer atom/residue information from the template.

    Coordinates, occupancy and B-factor come from the CREST/Open Babel PDB.
    Atom names, residue names, chain, element and charge come from template.

    Atom order is used for matching.
    """

    output = []
    atom_index = 0

    for line in model_lines:

        if not line.startswith(("ATOM  ", "HETATM")):
            output.append(line)
            continue

        atom_index += 1

        if atom_index not in template_atoms:
            raise ValueError(
                f"Conformer contains atom {atom_index}, "
                f"but template only contains "
                f"{len(template_atoms)} atoms."
            )

        template = template_atoms[atom_index]

        residue_number = residue_map[atom_index]
        residue_name = component_names[residue_number - 1]

        # --------------------------------------------------------
        # Validate element assignment
        # --------------------------------------------------------

        conformer_element = line[76:78].strip()
        template_element = template["element"].strip()

        if (
            conformer_element
            and template_element
            and conformer_element.upper() != template_element.upper()
        ):
            raise ValueError(
                f"Element mismatch at atom {atom_index}: "
                f"conformer={conformer_element}, "
                f"template={template_element}"
            )

        # --------------------------------------------------------
        # Construct PDB line
        # --------------------------------------------------------

        new_line = (
            line[:6] +
            line[6:12] +
            template["atom_name"] +
            template["altloc"] +
            f"{residue_name:>3}" +
            " " +                         # no chain
            f"{residue_number:>4}" +
            template["icode"] +
            line[27:76] +
            template["element"] +
            template["charge"] +
            "\n"
        )

        output.append(new_line)

    # ------------------------------------------------------------
    # Validate atom count
    # ------------------------------------------------------------

    if atom_index != len(template_atoms):
        raise ValueError(
            f"Atom count mismatch: conformer contains "
            f"{atom_index} atoms, template contains "
            f"{len(template_atoms)} atoms."
        )

    return output

def transfer_pdb_info(
    template_path,
    conformers_path,
    output_path,
):
    """
    Create a corrected multi-model PDB.

    Template:
        atom names + connectivity

    Conformers:
        coordinates + models

    Residue numbering is reconstructed from CONECT records.
    """

    # ------------------------------------------------------------
    # 1. Read template
    # ------------------------------------------------------------

    template_atoms, bonds = parse_template(template_path)

    print(f"Template atoms: {len(template_atoms)}")

    # ------------------------------------------------------------
    # 2. Determine molecular components
    # ------------------------------------------------------------

    components = find_connected_components(
        template_atoms,
        bonds,
    )

    print(f"Found {len(components)} connected components:")

    for i, component in enumerate(components, start=1):
        print(
            f"  Residue {i}: "
            f"{len(component)} atoms "
            f"(atoms {component[0]}–{component[-1]})"
        )

    # ------------------------------------------------------------
    # 3. Assign new residue numbers
    # ------------------------------------------------------------

    residue_map = make_residue_mapping(components)

    # ------------------------------------------------------------
    # 4. Determine residue names
    # ------------------------------------------------------------

    component_names = get_component_names(
        template_atoms,
        components,
    )

    print("\nResidues:")

    for i, name in enumerate(component_names, start=1):
        print(f"  {i}: {name}")

    # ------------------------------------------------------------
    # 5. Read conformer models
    # ------------------------------------------------------------

    models = parse_conformer_pdb(conformers_path)

    print(f"\nConformer models: {len(models)}")

    # ------------------------------------------------------------
    # 6. Update every model
    # ------------------------------------------------------------

    updated = []

    for model_number, model in enumerate(models, start=1):

        updated_model = update_model(
            model,
            template_atoms,
            residue_map,
            component_names,
        )

        updated.extend(updated_model)

    # ------------------------------------------------------------
    # 7. Write output
    # ------------------------------------------------------------

    with open(output_path, "w") as f:
        f.writelines(updated)

    print(f"\nWritten: {output_path}")

    return output_path