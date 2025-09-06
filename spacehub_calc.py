import pandas as pd
import numpy as np
import numbers

import warnings
warnings.filterwarnings("ignore")

class TwoBodyOrbit:
    
    npoints = 0
    df = None
    i = 0 #Primary mass object
    j = 1 #Secondary mass object

    #Setters and Update methods
    @classmethod
    def set_data(cls, value):
        cls.df = value

    @classmethod
    def set_ij(cls, i, j):
        cls.i = i
        cls.j = j

    @classmethod
    def set_npoints(cls):
        df = cls.df
        cls.numpoints = df["time"].nunique()

    def __init__(self, filename, i, j):
        self.df = load_spacehub_data(filename)
        self.i = i
        self.j = j
        self.set_npoints()



#-----Helper/Standalone Functions Below-----
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

def calc_h_vector(r, v):
    #get h vector (Murray-Dermott eq 2.129) at a single point in time
    x, y, z = r[0], r[1], r[2]
    vx, vy, vz = v[0], v[1], v[2]

    if isinstance(x, numbers.Number):
        hx = y*vz-z*vy
        hy = z*vx - x*vz
        hz = x*vy - y*vx

        return hx, hy, hz

    hx = []
    hy = []
    hz = []

    for t in range(0, len(x)):
        
        hx.append(y[t]*vz[t]-z[t]*vy[t])
        hy.append(z[t]*vx[t] - x[t]*vz[t])
        hz.append(x[t]*vy[t] - y[t]*vx[t])

    return hx, hy, hz

#-----------Getter/Modifier functions----------
#
#------------Helpers/Intermediates---------------



def get_h_vector(data, i, j):
    dx, dy, dz = distance(data, 'p', i, j)
    dvx, dvy, dvz = distance(data, 'v', i, j)
    hvec = calc_h_vector([dx, dy, dz], [dvx, dvy, dvz])

    return hvec

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
    #get center of mass
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

#---------Outputs/Keplerian Elements-----------

def get_L(data, i, j):
    mi = get_tot_mass(data, i)
    mj = get_tot_mass(data, j)
    dx, dy, dz = distance(data, 'p', i, j)
    dvx, dvy, dvz = distance(data, 'v', i, j)

    Lx, Ly, Lz = [], [], []
    for t in range(0, len(dx)):
        Lxi, Lyi, Lzi = calc_L(mi, mj, dx[t], dy[t], dz[t], dvx[t], dvy[t], dvz[t])
        Lx.append(Lxi)
        Ly.append(Lyi)
        Lz.append(Lzi)

    return Lx, Ly, Lz

def mag(vec):
    if isinstance(vec[0], numbers.Number):
        square = 0
        for a in vec:
            square += a**2

        return np.sqrt(square)

    mags = []
    for t in range(0, len(vec[0])):
        square = vec[0][t]**2 + vec[1][t]**2 + vec[2][t]**2
        

        mags.append(np.sqrt(square))

    return mags

#---Keplerian orbital elements---
# Size and Shape
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
    #Nearly Eq. 2.134 of Murray-Dermott
    mi = get_tot_mass(data, i)
    mj = get_tot_mass(data, j)
    dx, dy, dz = distance(data, 'p', i, j)
    dvx, dvy, dvz = distance(data, 'v', i, j)

    smas = []
    for t in range(0, len(dx)):
        a = calc_sma(mi + mj, dx[t], dy[t], dz[t], dvx[t], dvy[t], dvz[t])
        smas.append(a)

    return smas

def get_scalar_e(data, i, j):
    #Eq 2.135 of Murray-Dermott
    mi = get_tot_mass(data, i)
    mj = get_tot_mass(data, j)
    mu = (mi+mj) #* G.value
    dx, dy, dz = distance(data, 'p', i, j)
    dvx, dvy, dvz = distance(data, 'v', i, j)
    hvec = calc_h_vector([dx, dy, dz], [dvx, dvy, dvz])
    h = mag(hvec)
    a = get_sma(data, i, j)
    try:
        return np.sqrt(1-(h**2)/(mu * a))
    except TypeError:
        e = []
        for t in range(0, len(a)):
            e.append(np.sqrt(1-(h[t]**2)/(mu * a[t])))
        
        return e

def get_inclination(data, i, j):
    #Eq 2.134 of Murray-Dermott
    mi = get_tot_mass(data, i)
    mj = get_tot_mass(data, j)
    mu = (mi+mj) #* G.value
    dx, dy, dz = distance(data, 'p', i, j)
    dvx, dvy, dvz = distance(data, 'v', i, j)
    hvec = calc_h_vector([dx, dy, dz], [dvx, dvy, dvz])
    h = mag(hvec)
    try:
        return np.rad2deg(np.arccos(hvec[2]/h))
    except TypeError:
        I = []
        for t in range(0, len(h)):
            I.append(np.rad2deg(np.arccos(hvec[2][t]/h[t])))
        return I
    

def get_longitude_of_ascending_node(data, i, j):
    #requires list input, returns Omega, sinOmega, cosOmega
    mi = get_tot_mass(data, i)
    mj = get_tot_mass(data, j)
    mu = (mi+mj) #* G.value
    dx, dy, dz = distance(data, 'p', i, j)
    dvx, dvy, dvz = distance(data, 'v', i, j)
    hvec = calc_h_vector([dx, dy, dz], [dvx, dvy, dvz])
    h = mag(hvec)
    incl = get_inclination(data, i, j)

 

    sines = []
    cosines = []
    Omegas = []
    for t in range(0, len(h)):
        if hvec[2][t] > 0:
            hx = hvec[0][t]
            hy = -1 * hvec[1][t]
        else:
            hx = -1 * hvec[0][t]
            hy = hvec[1][t]

        sines.append(hx / (h[t]*np.sin(np.deg2rad(incl[t]))))
        cosines.append(hy / (h[t]*np.sin(np.deg2rad(incl[t]))))
        Omegas.append(np.rad2deg(np.arcsin(hx / (h[t]*np.sin(np.deg2rad(incl[t]))))))
    return Omegas, sines, cosines



def get_true_anomaly(data, i, j):
    a = get_sma(data, i, j)
    e = get_scalar_e(data, i, j)
    hvec = get_h_vector(data, i, j)
    h = mag(hvec)
    X, Y, Z = distance(data, 'p', i, j)
    R = mag([X, Y, Z])
    dX, dY, dZ = distance(data, 'v', i, j)

    #calculate rate of change of length of radius vector (Eq. 2.130)
    Rdot = []
    for t in range(0, len(a)):
        RdotRdot = X[t]*dX[t] +Y[t]*dY[t] + Z[t]* dZ[t] #2.128 
        V2 = dX[t]**2 + dY[t]**2 + dZ[t]**2

        Rdot.append(np.sign(RdotRdot) * np.sqrt(V2 - (h[t]**2 / R[t]**2)))


    sinf = []
    cosf = []
    f = []
    for t in range(0, len(a)):
        a1e2 = a[t] * (1-e[t]**2)
        sinft = (a1e2) * Rdot[t] / (h[t] * e[t])  
        cosft = (1/e[t]) * (  (a1e2/R[t]) - 1 )
        sinf.append(sinft)
        cosf.append(cosft)

        f.append(np.rad2deg(np.arcsin(sinft)))

    return f, sinf, cosf

def get_argument_of_periapsis(data, i, j):
    Omega, sinOmega, cosOmega = get_longitude_of_ascending_node(data, i, j)
    f, sinf, cosf = get_true_anomaly(data, i, j)
    incl = get_inclination(data, i, j)
    X, Y, Z = distance(data, 'p', i, j)
    R = mag([X,Y,Z])

    omega = []
    omega1 = []
    for t in range(0, len(f)):
        sinomegaft = Z[t] / (R[t]*np.sin(np.deg2rad(incl[t])))
        cosomegaft = (1/cosOmega[t]) * ((X[t]/R[t])+sinOmega[t]*sinomegaft*np.cos(np.deg2rad(incl[t])))
        
        fradianst = np.arcsin(sinf[t])

        omegat = np.rad2deg(np.arcsin(sinomegaft) - fradianst)
        omegatc = np.rad2deg(np.arccos(cosomegaft) - fradianst)
        
        omega.append(omegat)
        omega1.append(omegatc)

    return omega, omega1


#-------------Read/Write Operations--------------
# Load DefaultWriter output
def load_spacehub_data(filename):
    df = pd.read_csv(filename)
    add_norms(df)
    return df

