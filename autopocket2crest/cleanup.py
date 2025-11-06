import os
import shutil

def cleanup_temp_files(file_list):
    """Delete intermediate files generated during AutoPocket2CREST."""
    for f in file_list:
        if os.path.exists(f):
            try:
                os.remove(f)
                print(f"Removed: {f}")
            except Exception as e:
                print(f"Could not remove {f}: {e}")

def move_results_to_target(pdbid, target_dirs):
    """
    Move final results to shared target directories.
    target_dirs: dict like {"core": ".../core_set", "refined": ".../refined_set"}
    """
    source = os.path.join("data", "local", pdbid)
    moved = False
    for name, path in target_dirs.items():
        dest = os.path.join(path, pdbid)
        if os.path.exists(dest):
            for f in os.listdir(source):
                shutil.move(os.path.join(source, f), dest)
            print(f"Results moved to {name} set: {dest}")
            moved = True
            break

    if not moved:
        print(f"Warning: {pdbid} not found in any known target dataset. Files remain in {source}")

def safe_rmdir(path):
    """Remove a directory and its contents safely."""
    if os.path.exists(path) and os.path.isdir(path):
        shutil.rmtree(path)
        print(f"Removed directory: {path}")
