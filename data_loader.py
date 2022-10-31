import re
import os

digit = r'\d+'

def get_metadata(file):
    name, ext = file.split('.')
    if ext == "Chn":
        return None
    unprocessed = name.split('_')
    iteration = 1 if len(unprocessed) == 2 else int(re.search(digit, unprocessed[2]).group(0))
    return unprocessed[0], unprocessed[1], iteration

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