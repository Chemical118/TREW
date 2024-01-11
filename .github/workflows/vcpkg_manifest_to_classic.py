import sys
import json

json_path = sys.argv[1]

with open(json_path, 'r') as f:
    data = json.load(f)
    dep_list = [i["name"] for i in data['dependencies']]
    print(*dep_list, sep='\n')
