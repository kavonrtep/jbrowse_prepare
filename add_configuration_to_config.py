# script to add configuration to json file
import json
import sys
# get positional argument- json file
json_file = sys.argv[1]

# load json file
with open(json_file, 'r') as f:
    data = json.load(f)
# add configuration to json file
# "configuration": {"logoPath": {"uri": "./elixir_150x48.svg"}},
data["configuration"] = {"logoPath": {"uri": "./elixir_150x48.svg"}}

# write to json file
with open(json_file, 'w') as f:
    json.dump(data, f, indent=4)
