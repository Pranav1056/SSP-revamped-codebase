import numpy as np

def getMomentum(fileName):
    file = open(fileName).read().split()
    print(file)
    r = np.array([float(file[3]), float(file[4]), float(file[5])])
    v = np.array([float(file[6]), float(file[7]), float(file[8])])

    v = v*(365.2563835/(2*np.pi))

    return np.cross(r, v)

print(getMomentum("D:\Ji_OD\JiInput.txt"))