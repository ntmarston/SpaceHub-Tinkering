import pandas as pd
import numpy as np
import numbers
from matplotlib import pyplot as plt


import warnings
warnings.filterwarnings("ignore")

class TwoBodyOrbit:
    
    npoints = 0
    data = None
    i = 0 #Primary mass object
    j = 1 #Secondary mass object

    #Orbit State Vectors
    time = None
    R_vec = None
    V_vec = None
    magR = None
    hvec = None
    e_vec = None
    N_vec = None
    magN = None

    eccentricity = None
    semiMajorAxis = None
    inclination_rad = None
    inclination_deg = None


    LongitudeAscendingNode = None
    LongitudeAscendingNode_deg = None
    sinOmega = None
    cosOmega = None

    true_anomaly_rad = None
    true_anomaly_deg = None
    eccentric_anomaly_rad = None
    eccentric_anomaly_deg = None
    sinf = None
    cosf = None

    argument_of_periapsis_rad = None
    argument_of_periapsis_deg = None

    time_of_pericenter_passage = None
    #Setters and Update methods
    
    def set_data(self, value):
        self.data = value
    
    def set_time(self):
        col = self.data["time"]
        self.time = col.drop_duplicates().to_list()
    
    def set_ij(self, i, j):
        self.i = i
        self.j = j

    
    def set_npoints(self):
        df = self.data
        timesteps = df["time"].nunique()
        self.npoints = timesteps

    def set_R_and_V(self):
        data, i, j = self.data, self.i, self.j
        X, Y, Z = distance(data, 'p', i, j)
        magR = mag([X, Y, Z])
        self.R_vec = [X, Y, Z]
        self.magR = magR

        vx, vy, vz = distance(data, 'v', i, j)
        self.V_vec = [vx, vy, vz]

    def set_h_vector(self):
        i = self.i
        j = self.j
        data = self.data
        dx, dy, dz = distance(data, 'p', i, j)
        dvx, dvy, dvz = distance(data, 'v', i, j)
        self.hvec = calc_h_vector([dx, dy, dz], [dvx, dvy, dvz])

    
    def set_ecc_vector(self):
        i = self.i
        j = self.j
        data = self.data
        mi = get_tot_mass(data, i)
        mj = get_tot_mass(data, j)
        dx, dy, dz = distance(data, 'p', i, j)
        dvx, dvy, dvz = distance(data, 'v', i, j)

        ecx, ecy, ecz = [], [], []
        for t in range(0, len(dx)):
            ex, ey, ez = calc_ecc(mi + mj, dx[t], dy[t], dz[t], dvx[t], dvy[t], dvz[t])
            ecx.append(ex)
            ecy.append(ey)
            ecz.append(ez)
        
        self.e_vec = [ecx, ecy, ecz]
    
    def set_N_vector(self):
        i = self.i
        j = self.j
        data = self.data
        nxs, nys, nzs = [], [], []
        N = []
        for t in range(0, self.npoints):
            r = [self.R_vec[0][t], self.R_vec[1][t], self.R_vec[2][t]]
            v = [self.V_vec[0][t], self.V_vec[1][t], self.V_vec[2][t]]
            h = np.cross(r,v)
            khat = [0,0,1]
            kcrossh = np.cross(khat, h)
            nx, ny, nz = kcrossh[0], kcrossh[1], kcrossh[2]
            nxs.append(nx)
            nys.append(ny)
            nzs.append(nz)
            N.append(mag([nx, ny, nz]))

        self.N_vec = [nxs, nys, nzs]
        self.magN = N

    def set_sma(self):
        i = self.i
        j = self.j
        data = self.data
        #Nearly Eq. 2.134 of Murray-Dermott
        mi = get_tot_mass(data, i)
        mj = get_tot_mass(data, j)
        dx, dy, dz = distance(data, 'p', i, j)
        dvx, dvy, dvz = distance(data, 'v', i, j)

        smas = []
        for t in range(0, len(dx)):
            a = calc_sma(mi + mj, dx[t], dy[t], dz[t], dvx[t], dvy[t], dvz[t])
            smas.append(a)

        self.semiMajorAxis = smas

    
    def set_scalar_e(self):
        i = self.i
        j = self.j
        data = self.data
        #Eq 2.135 of Murray-Dermott
        mi = get_tot_mass(data, i)
        mj = get_tot_mass(data, j)
        mu = (mi+mj) #* G.value
        hvec = self.hvec
        h = mag(hvec)
        a = self.semiMajorAxis

        e = []

        for t in range(0, len(a)):
            e.append(np.sqrt(1-(h[t]**2)/(mu * a[t])))
        
        self.eccentricity = e

   
    def set_inclination(self):
        i = self.i
        j = self.j
        data = self.data
        #Eq 2.134 of Murray-Dermott
        mi = get_tot_mass(data, i)
        mj = get_tot_mass(data, j)
        mu = (mi+mj) #* G.value
        dx, dy, dz = distance(data, 'p', i, j)
        dvx, dvy, dvz = distance(data, 'v', i, j)
        hvec = calc_h_vector([dx, dy, dz], [dvx, dvy, dvz])
        h = mag(hvec)
        I = []
        for t in range(0, len(h)):
            I.append(np.arccos(hvec[2][t]/h[t]))
        
        self.inclination_rad = I
        self.inclination_deg = np.rad2deg(I)

    
    def set_longitude_of_ascending_node(self):
        i = self.i
        j = self.j
        data = self.data

        #returns Omega, sinOmega, cosOmega
        mi = get_tot_mass(data, i)
        mj = get_tot_mass(data, j)
        mu = (mi+mj) #* G.value
        hvec = self.hvec
        h = mag(hvec)
        incl = self.inclination_rad

    

        sines = []
        cosines = []
        Omegas = []
        for t in range(0, self.npoints):
            if hvec[2][t] > 0:
                hx = hvec[0][t]
                hy = -1 * hvec[1][t]
            else:
                hx = -1 * hvec[0][t]
                hy = hvec[1][t]

        
            
                
            sines.append(hx / (h[t]*np.sin(incl[t])))
            cosines.append(hy / (h[t]*np.sin(incl[t])))

            Omega_t = np.arcsin(hx / (h[t]*np.sin(incl[t])))
            if np.isnan(Omega_t):
                Omega_t = - np.arccos(hy / (h[t]*np.sin(incl[t]))) #negative sign is a hotfix
                
            if Omega_t < 0:
                Omega_t = 2*np.pi + Omega_t
            Omegas.append(Omega_t)
            


        #checkpoint
        for t in range(0, self.npoints):
            checksum = sines[t]**2 + cosines[t]**2
            assert np.abs(1-checksum) < 0.1, f"Checkpoint test failed in set_longitude_of_ascending_node. sin^2+cos^2 = {checksum}"
        
        self.sinOmega = sines
        self.cosOmega = cosines
        self.LongitudeAscendingNode = Omegas
        self.LongitudeAscendingNode_deg = [np.rad2deg(o) for o in Omegas]
    
    def set_true_anomaly(self):
        #Could not make the equations in Murray-Dermott cooperate with the way numpy trig functions work,
        #So I am using the eccentricity vector to calculate it
        evec = self.e_vec
        rvec = self.R_vec
        vvec = self.V_vec
        f = []
        fdeg = []
        cosf = []
        sinf = []
        for t in range(0, self.npoints):
            e = [evec[0][t], evec[1][t], evec[2][t]]
            r = [rvec[0][t], rvec[1][t], rvec[2][t]]
            v = [vvec[0][t], vvec[1][t], vvec[2][t]]
            edotr = scalar_product_xyz(e, r)
            mager = mag(e) * mag(r)
            rdotv = scalar_product_xyz(r, v)
            cosft = edotr/mager
            if cosft > 1: #Added because precision in spacehub sometimes yields cosf ~1+1e-7 which causes numpy to crash
                cosft = 1.0

            cosf.append(cosft)

            if rdotv > 0:
                ft = np.arccos(cosft)
            else:
                ft = 2*np.pi - np.arccos(cosft)
            
            
            sinf.append(np.sqrt(1-cosft**2))
            fdeg.append(np.rad2deg(ft))
            f.append(ft)

        self.true_anomaly_rad = f
        self.true_anomaly_deg = fdeg
        self.sinf = sinf
        self.cosf = cosf
    
    def set_argument_of_periapsis(self):
        evec = self.e_vec
        escalars = self.eccentricity
        peri = []
        perideg = []
        nvec = self.N_vec
        N = self.magN
        for t in range(0, self.npoints):
            Nv = [nvec[0][t], nvec[1][t], nvec[2][t]]
            e_vec = [evec[0][t], evec[1][t], evec[2][t]]
            e = escalars[t]
            Nt = N[t]
            if e_vec[2] >= 0:
                omega = np.arccos(np.dot(Nv, e_vec) / (Nt * e))
            else:
                omega = 2 * np.pi - np.arccos(np.dot(Nv, e_vec) / (Nt * e))
            peri.append(omega)
            perideg.append(np.rad2deg(omega))

        self.argument_of_periapsis_rad = peri
        self.argument_of_periapsis_deg = perideg

    def set_eccentric_anomaly(self):
        #Eq. 2.42 Murray-Dermott
        vvec = self.V_vec
        rvec = self.R_vec
        r_list = self.magR
        e_list = self.eccentricity
        a_list = self.semiMajorAxis
        eccentric_anomaly_rad = []
        eccentric_anomaly_deg = []
        for t in range(0, self.npoints):
            rv = [rvec[0][t], rvec[1][t], rvec[2][t]]
            vv = [vvec[0][t], vvec[1][t], vvec[2][t]]
            rdotv = scalar_product_xyz(rv, vv)
            r = r_list[t]
            e = e_list[t]
            a = a_list[t]
            arg = (1/e)*(-r/a+1)
            if rdotv > 0:
                E = np.arccos(arg)
            else:
                E = 2*np.pi - np.arccos(arg)
            eccentric_anomaly_rad.append(E)
            eccentric_anomaly_deg.append(np.rad2deg(E))
            

        self.eccentric_anomaly_rad = eccentric_anomaly_rad
        self.eccentric_anomaly_deg = eccentric_anomaly_deg

    def set_time_of_pericenter_passage(self):
        Gvalue = 1
        time = self.time
        eccentricity = self.eccentricity
        eccentric_anomaly = self.eccentric_anomaly_rad
        semiMajor = self.semiMajorAxis
        mi = get_tot_mass(self.data, self.i)
        mj = get_tot_mass(self.data, self.j)
        mu = (mi+mj) #* G.value
        taus = []
        for p in range(0, self.npoints):
            t = time[p]
            E = eccentric_anomaly[p]
            e = eccentricity[p]
            a = semiMajor[p]
            tau = t - (E - e*np.sin(E))/np.sqrt(mu * a**(-3))
            taus.append(tau)
        
        self.time_of_pericenter_passage = taus
    
    def __init__(self, filename, i, j):
        self.data = load_spacehub_data(filename)
        self.i = i
        self.j = j
        #set len and time
        self.set_npoints()
        self.set_time()
        #set R vec
        self.set_R_and_V()
        #set hvec
        self.set_h_vector()
        #set N vec
        self.set_N_vector()
        #solve eccentricity vector
        self.set_ecc_vector()
        #solve Semi-major axis
        self.set_sma()
        #solve eccentricity
        self.set_scalar_e()
        #solve inclination
        self.set_inclination()
        #solve longitude of ascending node
        self.set_longitude_of_ascending_node()
        #solve true anomaly
        self.set_true_anomaly()
        #solve argument of periapsis
        self.set_argument_of_periapsis()
        #solve eccentric anomaly
        self.set_eccentric_anomaly()
        #solve time of pericenter passage
        self.set_time_of_pericenter_passage()


    def __str__(self):
        outstr = ""
        outstr += (f"Initial (t=0) Conditions:" + "\n"
                   + f"a: {self.semiMajorAxis[0]}AU" + "\n"
                   + f"e: {self.eccentricity[0]}" + "\n"
                   + f"i: {self.inclination_deg[0]}deg" + "\n" 
                   + f"Longtiude of Ascending Node: {self.LongitudeAscendingNode_deg[0]}deg" + "\n"
                   + f"Argument of Periapsis: {self.argument_of_periapsis_deg[0]}deg" + "\n" 
                   + f"True Anomaly: {self.true_anomaly_deg[0]}deg" + "\n"
                   + f"Orbital Separation: {self.magR[0]}" + "\n")
        return outstr


    def plot_orbit_3panel(self):
        df, i, j = self.data, self.i, self.j
        fig = plt.figure(figsize=(15, 5))
        ax0 = fig.add_subplot(131, projection='3d')
        ax1 = fig.add_subplot(132)
        ax2 = fig.add_subplot(133)
        ax0.plot([0,1],[0],[0])
        ax0.plot([0],[0,1],[0])
        ax0.plot([0],[0],[0,1])
        ax0.scatter(df[df["id"]==j]["px"],df[df["id"]==j]["py"],df[df["id"]==j]["pz"], s=0.3, c='red', zorder=2)
        ax0.scatter(df[df["id"]==i]["px"],df[df["id"]==i]["py"],df[df["id"]==i]["pz"], s=1, c='blue', zorder=1)
        ax0.set_xlim(-8, 8)
        ax0.set_ylim(-8, 8)
        ax1.scatter(df[df["id"]==i]["px"], df[df["id"]==i]["py"], s=1, c='blue', zorder=1)
        ax1.scatter(df[df["id"]==j]["px"], df[df["id"]==j]["py"], s=0.3, c='red', zorder=2)
        ax1.grid(visible=True, zorder=-1)

        ax2.scatter(df[df["id"]==j]["px"], df[df["id"]==j]["pz"], s=0.3, c='red', zorder=2)
        ax2.scatter(df[df["id"]==i]["px"], df[df["id"]==i]["pz"], s=1, c='blue', zorder=1)
        ax2.grid(visible=True, zorder=-1)

        return fig, [ax0, ax1, ax2]
    
    def plot_keplerian_evolution_basic(self, xlim = []):
        #e, a, i, Omega, f, omega
        if len(xlim) < 2:
            xlim = [0, max(self.time)]

        fig, axes = plt.subplots(2,3, figsize=(15, 10))
        axs = axes.flatten()
        axs[0].plot(self.time, self.eccentricity)
        axs[0].set_ylim(-0.1, 1)
        axs[0].set_ylabel("Eccentricity")
        axs[0].set_xlim(xlim)

        axs[1].plot(self.time, self.inclination_deg)
        axs[1].set_xlim(xlim)
        axs[1].set_ylim(0, 360)
        axs[1].set_ylabel(r"Inclination $i$ (deg)")

        axs[2].plot(self.time, self.magR)
        axs[2].set_xlim(xlim)
        axs[2].set_title("R (AU)")
        axs[2].set_ylim(0, 10)
        axs[2].set_ylabel(r"Separation $||R||$ ($AU$)")

        axs[3].plot(self.time, self.LongitudeAscendingNode_deg)
        axs[3].set_title("")
        axs[3].set_xlim(xlim)
        axs[3].set_ylim(0, 360)
        axs[3].set_ylabel(r"Longitude of Ascending Node $\Omega$ (deg)")


        axs[4].plot(self.time, self.true_anomaly_deg)
        axs[4].set_xlim(xlim)
        axs[4].set_ylim(0, 360)
        axs[4].set_ylabel(r"True Anomaly $f$ (deg)")

        axs[5].plot(self.time, self.argument_of_periapsis_deg)
        axs[5].set_xlim(xlim)
        axs[5].set_ylim(0, 360)
        axs[5].set_ylabel(r"Argument of Periapsis $\omega$ (deg)")

        for ax in axs:
            ax.set_xlabel("$yr (2\pi)^{-1}$")
        return fig, axs

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

def get_tot_mass(data, tup):
    """Get the total mass of one or more particles (passed as tuple or int)

    Args:
        data (pd.DataFrame): The spacehub output dataframe created by load_spacehub_data
        param2 (int or tuple): The index/indices of the particles to consider

    Returns:
        float: The total mass of the particle(s)

    """
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
#-----------Getter/Modifier functions----------
#
#------------Helpers/Intermediates---------------

def get_h_vector(data, i, j):
    dx, dy, dz = distance(data, 'p', i, j)
    dvx, dvy, dvz = distance(data, 'v', i, j)
    hvec = calc_h_vector([dx, dy, dz], [dvx, dvy, dvz])

    return hvec


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

def scalar_product_xyz(vec1, vec2):
    #does not do lists
    return vec1[0]*vec2[0] + vec1[1]*vec2[1] + vec1[2]*vec2[2]

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





#-------------Read/Write Operations--------------
# Load DefaultWriter output
def load_spacehub_data(filename):
    #units: AU = 1, year = 2pi, G = 1
    df = pd.read_csv(filename)
    add_norms(df)
    return df

