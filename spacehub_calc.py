import pandas as pd
import numpy as np
import warnings
warnings.filterwarnings("ignore")

#Calc functions
def calc_norm(x, y, z):
    return np.sqrt(x ** 2 + y ** 2 + z ** 2)

def calc_ecc(m_tot, dx, dy, dz, dvx, dvy, dvz):
    u = m_tot * 1 #G.value
    v2 = dvx ** 2 + dvy ** 2 + dvz ** 2
    r = np.sqrt(dx ** 2 + dy ** 2 + dz ** 2)
    rv = dx * dvx + dy * dvy + dz * dvz
    ex = (dx * (v2 - u / r) - dvx * rv) / u
    ey = (dy * (v2 - u / r) - dvy * rv) / u
    ez = (dz * (v2 - u / r) - dvz * rv) / u

    return ex, ey, ez

def calc_sma(m_tot, dx, dy, dz, dvx, dvy, dvz):
    u = m_tot * 1 #G.value
    v2 = dvx ** 2 + dvy ** 2 + dvz ** 2
    r = np.sqrt(dx ** 2 + dy ** 2 + dz ** 2)
    return - u * r / (r * v2 - 2 * u)

def calc_angle(x1, y1, z1, x2, y2, z2):
    r1 = calc_norm(x1, y1, z1)
    r2 = calc_norm(x2, y2, z2)
    cos = (x1 * x2 + y1 * y2 + z1 * z2) / (r1 * r2)
    return np.arccos(cos)

def calc_L(m1, m2, dx, dy, dz, dvx, dvy, dvz):
    m_nu = m1 * m2 / (m1 + m2)
    Lx = dy * dvz - dz * dvy
    Ly = dz * dvx - dx * dvz
    Lz = dx * dvy - dy * dvx

    return m_nu * Lx, m_nu * Ly, m_nu * Lz

def distance(data, key, i, j):
    if type(i) is int:
        xi = data[data["id"]==i][key + 'x']
        yi = data[data["id"]==i][key + 'y']
        zi = data[data["id"]==i][key + 'z']
    elif type(i) is tuple:
        xi, yi, zi = get_com(data, key, i)
    else:
        print('wrong index type of i')

    if type(j) is int:
        xj = data[data["id"]==j][key + 'x']
        yj = data[data["id"]==j][key + 'y']
        zj = data[data["id"]==j][key + 'z']
    elif type(j) is tuple:
        xj, yj, zj = get_com(data, key, j)
    else:
        print('wrong index type of j')

    #de-index to avoid nans
    xdist, ydist, zdist = [], [], []
    for t in range(0, len(xi)):
        xdist.append(xi.iloc[t] - xj.iloc[t])
        ydist.append(yi.iloc[t] - yj.iloc[t])
        zdist.append(zi.iloc[t] - zj.iloc[t])
    return xdist, ydist, zdist

#Getter/Modifier functions
def get_tot_mass(data, tup):
    #get total mass of the system
    if type(tup) is int:
        return data[data["id"]==tup]["mass"][tup]
    else:
        mtot = 0

        for t in tup:
            mtot += data[data["id"]==t]["mass"][t]
        return mtot

def get_com(data, key, tup):
    mt = get_tot_mass(data, tup)

    x = 0
    y = 0
    z = 0

    for t in tup:
        x += data[data["id"]==t]['mass'] * data[data["id"]==t][key + 'x']
        y += data[data["id"]==t]['mass'] * data[data["id"]==t][key + 'y']
        z += data[data["id"]==t]['mass'] * data[data["id"]==t][key + 'z']

    return x/mt, y/mt, z/mt

def add_norms(data):
    #adds norm column to output, is inplace
    p_num = data["id"].nunique()
    for i in range(p_num):
        px = data['px']
        py = data['py']
        pz = data['pz']
        data['p'] = calc_norm(px, py, pz)

        vx = data['vx']
        vy = data['vy']
        vz = data['vz']
        data['v'] = calc_norm(vx, vy, vz)

def get_ecc(data, i, j):
    mi = get_tot_mass(data, i)
    mj = get_tot_mass(data, j)
    dx, dy, dz = distance(data, 'p', i, j)
    dvx, dvy, dvz = distance(data, 'v', i, j)

    ecc = []
    for t in range(0, len(dx)):
        e = calc_ecc(mi + mj, dx[t], dy[t], dz[t], dvx[t], dvy[t], dvz[t])
        ecc.append(e)
    return ecc

def get_sma(data, i, j):
    mi = get_tot_mass(data, i)
    mj = get_tot_mass(data, j)
    dx, dy, dz = distance(data, 'p', i, j)
    dvx, dvy, dvz = distance(data, 'v', i, j)

    smas = []
    for t in range(0, len(dx)):
        a = calc_sma(mi + mj, dx[t], dy[t], dz[t], dvx[t], dvy[t], dvz[t])
        smas.append(a)

    return smas

def get_L(data, i, j):
    mi = get_tot_mass(data, i)
    mj = get_tot_mass(data, j)
    dx, dy, dz = distance(data, 'p', i, j)
    dvx, dvy, dvz = distance(data, 'v', i, j)

    Ls = []
    for t in range(0, len(dx)):
        L = calc_L(mi, mj, dx[t], dy[t], dz[t], dvx[t], dvy[t], dvz[t])
        Ls.append(L)

    return Ls

def load_spacehub_data(filename):
    df = pd.read_csv("tutorial/hierarchical.txt")
    add_norms(df)
    return df