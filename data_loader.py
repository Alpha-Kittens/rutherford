import re

digit = r'\d+'

def get_metadata(file):
    name, ext = file.split(".")
    if ext == "Chn":
        return None
    unprocessed = name.split("_")
    iteration = 1 if len(unprocessed) == 2 else re.search(digit, unprocessed[2]).group(0)
    return unprocessed[0], unprocessed[1], iteration

def read_data(file, return_hist = False):
    
    metadata = get_metadata(file)
    if metadata is None:
        return None

    
