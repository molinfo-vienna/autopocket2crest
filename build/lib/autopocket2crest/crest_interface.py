import os
import subprocess


def compress_ranges(indices):
    """Compress integer list to compact range strings."""
    if not indices:
        return []

    indices = sorted(set(indices))

    ranges = []
    start = prev = indices[0]

    for n in indices[1:]:
        if n == prev + 1:
            prev = n
        else:
            ranges.append(
                f"{start}-{prev}" if start != prev else f"{start}"
            )
            start = prev = n

    ranges.append(
        f"{start}-{prev}" if start != prev else f"{start}"
    )

    return ranges


def generate_constraints(xyz_file, backbone_indices):
    """
    Generate constraints.inp file for CREST based on atom indices.

    Fails immediately if CREST cannot generate the constraint file.
    """

    constrain_str = ",".join(
        compress_ranges(backbone_indices)
    )

    cmd = [
        "crest",
        xyz_file,
        "--constrain",
        constrain_str,
    ]

    print("Generating CREST constraints...")
    print("Executing:", " ".join(cmd))

    result = subprocess.run(
        cmd,
        capture_output=True,
        text=True,
    )

    if result.returncode != 0:
        raise RuntimeError(
            "CREST failed while generating constraints.\n\n"
            f"STDOUT:\n{result.stdout}\n\n"
            f"STDERR:\n{result.stderr}"
        )

    if not os.path.exists(".xcontrol.sample"):
        raise RuntimeError(
            "CREST completed, but .xcontrol.sample was not created. "
            "Cannot generate constraints.inp."
        )

    os.replace(".xcontrol.sample", "constraints.inp")

    print("Constraint file created.")

    return "constraints.inp"


def run_crest(
    xyz_file,
    constraint_file=None,
    charge=0,
    temp="310",
    lvl_of_theory="gfnff",
    extra_crest_args="-squick",
):
    """
    Run CREST conformer search.

    Raises RuntimeError if CREST exits unsuccessfully.
    """

    cmd = [
        "crest",
        xyz_file,
        f"-{lvl_of_theory}",
        "-chrg",
        str(charge),
        "-gbsa",
        "h2o",
    ]

    # Handle extra arguments such as "-squick"
    if extra_crest_args:
        cmd.extend(extra_crest_args.split())

    cmd.extend([
        "--temp",
        str(temp),
        "--legacy",
    ])

    if constraint_file and os.path.exists(constraint_file):
        cmd.extend([
            "-cinp",
            constraint_file,
        ])

    print("Executing:", " ".join(cmd))

    with open("crest.out", "w") as stdout_file:

        result = subprocess.run(
            cmd,
            stdout=stdout_file,
            stderr=subprocess.STDOUT,
            text=True,
        )

    if result.returncode != 0:
        raise RuntimeError(
            f"CREST failed with exit code {result.returncode}. "
            "See crest.out for details."
        )

    # Don't claim success unless CREST actually produced output.
    if not os.path.exists("crest_conformers.xyz"):
        raise RuntimeError(
            "CREST returned successfully, but "
            "crest_conformers.xyz was not produced. "
            "See crest.out for details."
        )

    print("CREST conformer search completed successfully.")

    return "crest.out"