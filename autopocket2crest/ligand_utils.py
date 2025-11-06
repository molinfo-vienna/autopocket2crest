def get_ligand_name(mol2_file, pdb_file):
    """Detect ligand residue name using MOL2 and fallback PDB scanning."""
    typical_cofactor_names = {
        "FAD","FMN","NAD","NADP","NADH","NADPH","ATP","ADP","GTP","GDP","CoA","CoQ",
        "HEM","HEME","PQQ","TPP","PLP","THF","FOL","VB6","VB12","VB1","VB2",
        "HOH","H2O","WAT","SO4","SO3","SO2","PO4","PO3","PO2","PO","CL","BR","I",
        "MG","CA","ZN","FE","CU","MN","NI","CO","MO","CR","V","W","SE","AS","AL",
        "BA","CD","HG","PB","PT","AU","AG","LI","NA","K","RB","CS","BE","SR","ZR",
        "Y","LA","CE","PR","ND","PM","SM","EU","GD","TB","DY","HO","ER","TM","YB",
        "LU","HF","TA","RE","OS","IR","BI","TH","PA","U","NP","PU","AM","CM","BK",
        "CF","ES","FM","MD","NO","LR","RF","DB","SG","BH","HS","MT","DS","RG","CN",
        "NH","FL","MC","LV","TS","OG"
    }

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
        with open(pdb_file, 'r') as pdb:
            for line in pdb:
                if line.startswith("HETATM"):
                    res = line[17:20].strip()
                    if res not in typical_cofactor_names:
                        ligand_name = res
                        break

    if not ligand_name:
        print(f"Warning: No ligand name found for {mol2_file}, using 'UNNAMED'")
        ligand_name = "UNNAMED"

    return ligand_name
