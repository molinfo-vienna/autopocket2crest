import os

def compress_ranges(indices):
    """Compress integer list to compact range strings."""
    if not indices: return []
    indices = sorted(set(indices))
    ranges, start, prev = [], indices[0], indices[0]
    for n in indices[1:]:
        if n == prev + 1:
            prev = n
        else:
            ranges.append(f"{start}-{prev}" if start != prev else f"{start}")
            start = prev = n
    ranges.append(f"{start}-{prev}" if start != prev else f"{start}")
    return ranges

def generate_constraints(xyz_file, backbone_indices):
    """Generate constraints.inp file for CREST based on atom indices."""
    constrain_str = ",".join(compress_ranges(backbone_indices))
    os.system(f"crest {xyz_file} --constrain {constrain_str}")
    if os.path.exists(".xcontrol.sample"):
        os.rename(".xcontrol.sample", "constraints.inp")
    print("Constraint file created.")
    return "constraints.inp"

def run_crest(xyz_file, constraint_file=None, charge=0):
    """Run CREST conformer search."""
    cmd = f"crest {xyz_file} -gfnff -chrg {charge} -gbsa h2o -squick --temp 310 --legacy"
    if constraint_file and os.path.exists(constraint_file):
        cmd += f" -cinp {constraint_file}"
    cmd += " > crest.out"
    print("Executing:", cmd)
    os.system(cmd)
    return "crest.out"
