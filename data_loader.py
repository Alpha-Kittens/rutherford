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
    #print (file)
    name, ext = file.split('/')[-1].split('.')
    if ext != "Spe":
        return None
    unprocessed = name.split('_')
    if len(unprocessed) == 3 and unprocessed[2] == 'bad':
        return None
    iteration = 1 if len(unprocessed) == 2 else int(re.search(digit, unprocessed[2]).group(0))
    return str.lower(unprocessed[0]), int(unprocessed[1]), iteration

def get_file_info(fp):
    """
    Given file path pointing to a .Spe file, returns total time and and a histogram of channel counts as a 2-tuple.
    """
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


def read_data(file, metadata = None):
    """
    Given a file path poitning to a .Spe file, returns a dictionary `data` with relevant informaiton.
    Keys:
        * `target`: foil in use. `empty`, `gold`, `2gold`, or `iron`. 
        * `angle`: angle of detector
        * `iteration`: iteration number. 
        * `time`: total time of scan, truncated to seconds. 
        * `counts`: histogram of channel data
        * `cps`: counts per second
    """
    data = {}
    metadata = get_metadata(file) if not metadata else metadata
    if metadata is None:
        return None
    data['target'], data['angle'], data['iteration'] = metadata
    data['time'], data['histogram'] = get_file_info(file)
    data['counts'] = sum(data['histogram'])
    data['cps'] = data['counts'] / data['time']

    return data

def add_data(data, file, metadata = None):
    """
    Adds data associated with a file to a dictionary.
    Arguments:
        * `data`: Dictionary where info is stored
        * `file`: File path pointing to a .Spe file
    Returns:
        * Nothing. Rather, adds data to dictionary as follows:
            -Key: metadata of `file`. i.e., foil, angle, iteration. 
            -Value: `time`, `histogram` in file. 
    """
    entry = read_data(file, metadata = metadata)
    key = (entry['target'], entry['angle'], entry['iteration'])
    if key in data.keys():
        key = (entry['target'], entry['angle'], entry['iteration'] + 1)
    data[key] = entry['time'], entry['histogram']

def any_in(list, inlist):
    for entry in list:
        if entry in inlist:
            return True
    return False
def multiple_in(list, inlist):
    """
    Returns `True` if if all values of `list` (first argument) are present in `inlist` (second argument), `False` otherwise. 
    """
    success = True
    for entry in list:
        success = success and entry in inlist
    return success

def recursive_read(data, folder, require = [], reject = [], condition = lambda x : True):
    """
    Reads all entries in a given folder into `data` dictionary. Can specify metadata details to require or reject, as well as more specific inclusion conditions. 
    Arguments:
        * `data` (dict): dictionary to which data entries will be added. Keys and values will be same as specified in `add_data`.  
        * `folder` (string): file path of folder to read through
        * `require` (array): Any entries in here must be present in metadata of a file for the file to be read into the dictionary. Default: []
        * `reject` (array): If any entries in here are present in metadata of file, file will not be read into the dictionary. Default: []
        * `condition` (function): If specified, `condition(metadata)` must evaluate to `True` in order for file to be read into dictionary. Default: always True
    Returns:
        * Nothing. Rather, adds data to dictionary as follows:
            -Keys: metadata of all files read into dictionary. i.e., foil, angle, iteration. 
            -Values: `time`, `histogram` associated with file.    
    """
    for entry in os.listdir(folder):
        path = str(folder) + "/" + str(entry)
        if not os.path.isdir(path):
            if "unknown" not in entry and "sus" not in entry: 
                #print(path)
                metadata = get_metadata(entry)                    
                #print (metadata)
                if metadata is not None and multiple_in(require, metadata[0:-1]) and (reject == [] or not any_in(reject, metadata[0:-1])) and condition(metadata):
                    if metadata in data.keys():
                    
                        i = metadata[2]
                        while (metadata[0], metadata[1], i) in data.keys():
                            i += 1
                        print ("Warning: multiple files with the same metadata were found. Newer entry will be added with an increased iteration.\nMetadata: "+str(metadata)+" -> "+str((metadata[0], metadata[1], i)))
                        metadata = (metadata[0], metadata[1], i)
                    add_data(data, path, metadata = metadata)
        else:
            #print (entry)
            #print (folder)
            #if "Calibration" not in path:
            recursive_read(data, path, require = require, reject = reject, condition = condition) # for some reason, os.path.join isn't doing what I expect

import numpy as np

def iterationless(data):
    """
    Given `data` dict as returned by `add_data` and `recursive_read`, returns new dictionary where all 
    values which have the same `foil` and `angle` metadata are combined into an array under new key `(foil, angle)`. 
    """
    keys = list(data.keys())[:]
    new_data = {}
    for (foil, angle, iter) in keys:
        if (foil, angle) not in new_data.keys():
            new_data[(foil, angle)] = []
        #else:
            #print ("added new one")
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
        
