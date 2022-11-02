import re
import os

digit = r'\-?\d+'

element_map = {
    'gold' : 79,
    'iron' : 26,
    'titanium' : 22, 
}

def Z(element):
    return element_map[element]

def get_metadata(file):
    name, ext = file.split('/')[-1].split('.')
    if ext == "Chn":
        return None
    unprocessed = name.split('_')
    iteration = 1 if len(unprocessed) == 2 else int(re.search(digit, unprocessed[2]).group(0))
    return str.lower(unprocessed[0]), int(unprocessed[1]), iteration

def get_file_info(fp):
    with open(fp) as file:
        lines = file.readlines()
    times = lines[9].split(' ')
    time = int(times[0])
    print(time)
    counts = []
    for i in range(12, 2060):
        digit = r'\d+'

        count = int(re.search(digit, lines[i]).group(0))

        counts.append(count)

    return time, counts


def read_data(file):
    
    data = {}
    metadata = get_metadata(file)
    if metadata is None:
        return None
    data['target'], data['angle'], data['iteration'] = metadata
    data['time'], data['histogram'] = get_file_info(file)
    data['counts'] = sum(data['histogram'])
    data['cps'] = data['counts'] / data['time']

    return data

def add_data(data, file):

    entry = read_data(file)
    data[(entry['target'], entry['angle'], entry['iteration'])] = entry['time'], entry['histogram']

def recursive_read(data, folder):
    
    for entry in os.listdir(folder):
        if os.path.isfile(os.path.join(entry, folder)):
            add_data(data, entry)
        else:
            recursive_read(data, entry)

def iterationless(data):
    for key in data:
        if (key[0], key[1]) not in data.keys():
            data[(key[0], key[1])] = 0
        data[(key[0], key[1])] = data[(key[0], key[1])][0] + data[key][0], data[(key[0], key[1])][1] + data[key][1], 

        
