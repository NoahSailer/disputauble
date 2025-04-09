import yaml sys
with open(f'{str(sys.argv[1])}.yaml', 'r') as file: addData = yaml.safe_load(file)
with open(f'{str(sys.argv[2])}.yaml', 'r') as file: oldData = yaml.safe_load(file)
addData['output'] = oldData['output']
with open(f'{str(sys.argv[1])}_{str(sys.argv[2])}.yaml', 'w') as file: yaml.dump(yamlData,file)