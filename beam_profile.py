import os
from data_loader import read_data
import matplotlib.pyplot as plt
import re

folder = 'beam_profile/'

files = os.listdir(folder)

angles = []
cps = []
for file_name in files:
    fp = folder + file_name

    information = read_data(fp)

    angles.append(information['angle'])
    cps.append(information['cps'])

plt.scatter(angles, cps, marker='o')
plt.show()


fp = folder + 'empty_0.Spe'

'''
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

time, counts = get_file_info(fp)

plt.plot(counts)
plt.show()
'''