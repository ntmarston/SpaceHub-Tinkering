import pandas as pd
import numpy as np
import matplotlib as plt
const_G = 1
const_yr = 2 * np.pi
const_Ms = 1
const_deg = np.pi / 180
const_Myr = const_yr * 1e6
const_Gyr = const_yr * 1e9
const_day = const_yr / 365.25636042
const_hr = const_day / 24
const_sec = const_hr / 3600

const_au = 1
const_km = const_au / 149597870.7
const_pc = const_au * 648000 / np.pi
const_Rs = const_km * 6.957e5
const_Rj = const_km * 69911
const_kms = const_km / const_sec
const_Me = const_Ms * 3.003e-6
const_Mj = const_Me * 317.8
const_C = const_kms * 299792.458


def calc_norm(x, y, z):
    return np.sqrt(x ** 2 + y ** 2 + z ** 2)


def get_norm(p_num, data):
    for i in range(p_num):
        px = data['px' + str(i)]
        py = data['py' + str(i)]
        pz = data['pz' + str(i)]
        data['p' + str(i)] = calc_norm(px, py, pz)

        vx = data['vx' + str(i)]
        vy = data['vy' + str(i)]
        vz = data['vz' + str(i)]
        data['v' + str(i)] = calc_norm(vx, vy, vz)


def load_spacehub_data(fname):
    data = pd.read_csv(fname)
    dic = {name: data[name].to_numpy() for name in data.columns}
    particle_num = len(np.unique(dic['id']))

    splited_dic = {}
    splited_dic['time'] = dic['time'][dic['id'] == 0]
    for i in range(particle_num):
        for key in dic.keys():
            if key != 'time' and key != 'id':
                splited_dic[key + str(i)] = dic[key][dic['id'] == i]

    get_norm(particle_num, splited_dic)

    return particle_num, splited_dic


def calc_ecc(m_tot, dx, dy, dz, dvx, dvy, dvz):
    u = m_tot * const_G
    v2 = dvx ** 2 + dvy ** 2 + dvz ** 2
    r = np.sqrt(dx ** 2 + dy ** 2 + dz ** 2)
    rv = dx * dvx + dy * dvy + dz * dvz
    ex = (dx * (v2 - u / r) - dvx * rv) / u
    ey = (dy * (v2 - u / r) - dvy * rv) / u
    ez = (dz * (v2 - u / r) - dvz * rv) / u
    return ex, ey, ez


def calc_L(m1, m2, dx, dy, dz, dvx, dvy, dvz):
    m_nu = m1 * m2 / (m1 + m2)
    Lx = dy * dvz - dz * dvy
    Ly = dz * dvx - dx * dvz
    Lz = dx * dvy - dy * dvx

    return m_nu * Lx, m_nu * Ly, m_nu * Lz


def calc_sma(m_tot, dx, dy, dz, dvx, dvy, dvz):
    u = m_tot * const_G
    v2 = dvx ** 2 + dvy ** 2 + dvz ** 2
    r = np.sqrt(dx ** 2 + dy ** 2 + dz ** 2)
    return - u * r / (r * v2 - 2 * u)


def get_tot_mass(data, tup):
    if type(tup) is int:
        return data['mass'+str(tup)]
    else:
        mtot = 0

        for t in tup:
            mtot += data['mass' + str(t)]
        return mtot


def get_com(data, key, tup):
    mt = get_tot_mass(data, tup)

    x = 0
    y = 0
    z = 0

    for t in tup:
        x += data['mass' + str(t)] * data[key + 'x' + str(t)]
        y += data['mass' + str(t)] * data[key + 'y' + str(t)]
        z += data['mass' + str(t)] * data[key + 'z' + str(t)]

    return x/mt, y/mt, z/mt


def distance(data, key, i, j):
    if type(i) is int:
        xi = data[key + 'x' + str(i)]
        yi = data[key + 'y' + str(i)]
        zi = data[key + 'z' + str(i)]
    elif type(i) is tuple:
        xi, yi, zi = get_com(data, key, i)
    else:
        print('wrong index type of i')

    if type(j) is int:
        xj = data[key + 'x' + str(j)]
        yj = data[key + 'y' + str(j)]
        zj = data[key + 'z' + str(j)]
    elif type(j) is tuple:
        xj, yj, zj = get_com(data, key, j)
    else:
        print('wrong index type of j')

    return xi - xj, yi - yj, zi - zj


def get_ecc(data, i, j):
    mi = get_tot_mass(data, i)
    mj = get_tot_mass(data, j)
    dx, dy, dz = distance(data, 'p', i, j)
    dvx, dvy, dvz = distance(data, 'v', i, j)

    return calc_ecc(mi + mj, dx, dy, dz, dvx, dvy, dvz)


def get_sma(data, i, j):
    mi = get_tot_mass(data, i)
    mj = get_tot_mass(data, j)
    dx, dy, dz = distance(data, 'p', i, j)
    dvx, dvy, dvz = distance(data, 'v', i, j)

    return calc_sma(mi + mj, dx, dy, dz, dvx, dvy, dvz)


def get_L(data, i, j):
    mi = get_tot_mass(data, i)
    mj = get_tot_mass(data, j)
    dx, dy, dz = distance(data, 'p', i, j)
    dvx, dvy, dvz = distance(data, 'v', i, j)
    return calc_L(mi, mj, dx, dy, dz, dvx, dvy, dvz)


def calc_angle(x1, y1, z1, x2, y2, z2):
    r1 = calc_norm(x1, y1, z1)
    r2 = calc_norm(x2, y2, z2)
    cos = (x1 * x2 + y1 * y2 + z1 * z2) / (r1 * r2)
    return np.arccos(cos)