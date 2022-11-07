import re
import os

digit = r'\-?\d+'

element_map = {
    'gold' : 79,
    'iron' : 26,
    'titanium' : 22,
    '2gold' : 79,
    '3gold' : 79,
}

def Z(element):
    return element_map[element]

def get_metadata(file):
    print (file)
    name, ext = file.split('/')[-1].split('.')
    if ext != "Spe":
        return None
    unprocessed = name.split('_')
    if len(unprocessed) == 3 and unprocessed[2] == 'bad':
        return None
    iteration = 1 if len(unprocessed) == 2 else int(re.search(digit, unprocessed[2]).group(0))
    return str.lower(unprocessed[0]), int(unprocessed[1]), iteration

def get_file_info(fp):
    with open(fp) as file:
        lines = file.readlines()
    times = lines[9].split(' ')
    time = int(times[0])
    #print(time)
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
    key = (entry['target'], entry['angle'], entry['iteration'])
    if key in data.keys():
        key = (entry['target'], entry['angle'], entry['iteration'] + 1)
    data[key] = entry['time'], entry['histogram']

def multiple_in(list, inlist):
    #print (list)
    #print (inlist)
    success = True
    for entry in list:
        success = success and entry in inlist
    #print (success)
    return success

def recursive_read(data, folder, require = [], reject = [], condition = lambda x : True):
    #print (os.listdir(folder))
    for entry in os.listdir(folder):
        #print (path := str(folder) + "/" + str(entry))
        path = str(folder) + "/" + str(entry)
        if not os.path.isdir(path):
            if "unknown" not in entry and "sus" not in entry: 
                metadata = get_metadata(entry)
                #print (metadata)
                if metadata is not None and multiple_in(require, metadata[0:-1]) and (reject == [] or not multiple_in(reject, metadata[0:-1])) and condition(metadata):
                    add_data(data, path)
        else:
            #print (entry)
            #print (folder)
            recursive_read(data, path, require = require, reject = reject) # for some reason, os.path.join isn't doing what I expect

import numpy as np

def iterationless(data):
    keys = list(data.keys())[:]
    new_data = {}
    for (foil, angle, iter) in keys:
        if (foil, angle) not in data.keys():
            new_data[(foil, angle)] = []
        new_data[(foil, angle)].append(tuple(data[(foil, angle, iter)]))
    return new_data
        #print (data[(foil, angle)])
    #print (data.keys())

""" old version--decided cannot average between histograms
def iterationless(data):
    for key in data:
        if (key[0], key[1]) not in data.keys():
            data[(key[0], key[1])] = 0
        data[(key[0], key[1])] = data[(key[0], key[1])][0] + data[key][0], data[(key[0], key[1])][1] + data[key][1]
"""
        
