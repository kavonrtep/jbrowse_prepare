#!/usr/bin/env python3

# repeats.csv contain genome name and then path to directory.
# set mask using setfact to rwx recursively on all directories

import os

with open("repeats.csv") as f:
    for line in f:
        # line could be empty
        try:
            genome, path = line.strip().split("\t")
        except ValueError:
            continue
        print(f"Setting mask on {path}")
        # firt check what is the mask setting on the directory
        cmd = f"getfacl {path}"
        # run the command and capture the output
        mask = os.popen(cmd).read()
        # extract mask line - it starts with `mask`
        mask = [l for l in mask.split("\n") if l.startswith("mask")]
        if len(mask) == 0:
            print("No mask found")
            continue
        mask = mask[0].split(":")[-1]
        print(f"Mask is {mask}")
        if mask == "rwx":
            print("Mask is already set")
            continue
        cmd = f"setfacl -R -m mask::rwx {path}"
        print(cmd)
        os.system(cmd)
        print(f"Mask set on {path}")
        print("")
