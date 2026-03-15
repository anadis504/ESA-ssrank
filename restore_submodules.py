#!/usr/bin/env python3
import os
import subprocess
import configparser

GITMODULES_FILE = ".gitmodules"
HASHES_FILE = "submodule_hashes.txt"

# ---- parse module_hashes ----
hashes = {}
with open(HASHES_FILE) as f:
    for line in f:
        if not line.strip():
            continue
        parts = line.strip().split()
        path = parts[0]
        commit_hash = parts[-1]
        hashes[path] = commit_hash

# ---- parse .gitmodules ----
config = configparser.ConfigParser()
config.read(GITMODULES_FILE)

submodules = {}
for section in config.sections():
    path = config[section]["path"]
    url = config[section]["url"]
    submodules[path] = url

# ---- helper function to clone and checkout a repo ----
def clone_and_checkout(url, path, commit_hash):
    if not os.path.exists(path):
        print(f"Cloning {url} into {path}...")
        subprocess.run(["git", "clone", url, path], check=True)
    os.chdir(path)
    subprocess.run(["git", "fetch", "--all"], check=True)
    subprocess.run(["git", "checkout", commit_hash], check=True)
    # update nested submodules if any
    subprocess.run(["git", "submodule", "update", "--init", "--recursive"], check=True)
    os.chdir("..")

# ---- main loop ----
for path, url in submodules.items():
    print(f"Processing submodule {path}")
    if path not in hashes:
        print(f"Warning: no hash for {path} in {HASHES_FILE}, skipping...")
        continue
    commit_hash = hashes[path]
    clone_and_checkout(url, path, commit_hash)

print("All submodules restored at the specified commits!")
