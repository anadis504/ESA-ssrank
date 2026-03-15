#!/bin/bash
# This script restores all submodules to exact commits.

if [ ! -f submodule_hashes.txt ]; then
  echo "Error: submodule_hashes.txt not found!"
  exit 1
fi

while read path url hash; do
  echo "Restoring $path at $hash..."
  mkdir -p "$path"
  cd "$path"
  git init
  git remote add origin "$url"
  git fetch --depth 1 origin "$hash"
  git checkout "$hash"
  cd - > /dev/null
done < submodule_hashes.txt

echo "All submodules restored!"