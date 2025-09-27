import pandas as pd
import numpy as np
import numbers
from matplotlib import pyplot as plt
from matplotlib.ticker import FormatStrFormatter
import matplotlib.animation as animation
from astropy import units as u
from astropy.constants import G, c
import matplotlib as mpl
mpl.rcParams['animation.embed_limit'] = 250


import warnings
warnings.filterwarnings("ignore")

class TwoBodyOrbit:
    
    npoints = 0
    data = None
    i = 0 #Primary mass object
    j = 1 #Secondary mass object

    M_i = None
    M_j = None


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

    c0 = None
    points_per_orbit = None
    points_per_orbit_avg = None

    #Setters and Update methods
    
    def set_data(self, value):
        self.data = value.dropna()
    
    def set_npoints(self):
        timesteps0 = len(self.data[self.data["id"]==self.i]["time"])
        timesteps1 = len(self.data[self.data["id"]==self.j]["time"])
        self.npoints = min(timesteps0, timesteps1)

    def set_time(self):
        timesteps_i = len(self.data[self.data["id"]==self.i]["time"])
        timesteps_j = len(self.data[self.data["id"]==self.j]["time"])
        if timesteps_i < timesteps_j:
            col = self.data[self.data["id"]==self.i]["time"]
        else: #If they have the same number of timesteps, or if j has more
            col = self.data[self.data["id"]==self.j]["time"] 

        
        self.time = col.to_list()
        
    def set_ij(self, i, j):
        self.i = i
        self.j = j

    def set_masses(self):
        mi = get_tot_mass(self.data, self.i)
        mj = get_tot_mass(self.data, self.j)
        self.M_i = mi
        self.M_j = mj

    def set_R_and_V(self):
        data, i, j = self.data, self.i, self.j
        X, Y, Z = distance(data, 'p', i, j, self.npoints)
        magR = mag([X, Y, Z])
        self.R_vec = [X, Y, Z]
        self.magR = magR

        vx, vy, vz = distance(data, 'v', i, j, self.npoints)
        self.V_vec = [vx, vy, vz]

    def set_h_vector(self):
        i = self.i
        j = self.j
        data = self.data
        dx, dy, dz = distance(data, 'p', i, j, self.npoints)
        dvx, dvy, dvz = distance(data, 'v', i, j, self.npoints)
        self.hvec = calc_h_vector([dx, dy, dz], [dvx, dvy, dvz])
   
    def set_ecc_vector(self):
        i = self.i
        j = self.j
        data = self.data
        mi = get_tot_mass(data, i)
        mj = get_tot_mass(data, j)
        dx, dy, dz = distance(data, 'p', i, j, self.npoints)
        dvx, dvy, dvz = distance(data, 'v', i, j, self.npoints)

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
                try:
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
                except IndexError:
                    print(f"{t}, {self.npoints}, {len(self.time)}, {len(self.R_vec[0])}, {len(self.R_vec[1])}, {len(self.R_vec[2])}")
            self.N_vec = [nxs, nys, nzs]
            self.magN = N

    rdebug = []
    v2debug = []
    def set_sma(self):
        i = self.i
        j = self.j
        data = self.data
        #Nearly Eq. 2.134 of Murray-Dermott
        mi = get_tot_mass(data, i)
        mj = get_tot_mass(data, j)
        dx, dy, dz = distance(data, 'p', i, j, self.npoints)
        dvx, dvy, dvz = distance(data, 'v', i, j, self.npoints)
        
        smas = []
        for t in range(0, len(dx)):
            a = calc_sma(mi + mj, dx[t], dy[t], dz[t], dvx[t], dvy[t], dvz[t])
            smas.append(a)
            self.rdebug.append(np.sqrt(dx[t] ** 2 + dy[t] ** 2 + dz[t] ** 2))
            self.v2debug.append(dvx[t] ** 2 + dvy[t] ** 2 + dvz[t] ** 2)
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
        dx, dy, dz = distance(data, 'p', i, j, self.npoints)
        dvx, dvy, dvz = distance(data, 'v', i, j, self.npoints)
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
    
    def set_c0(self, deviation_check = True): #find the constant c0 defined in Peters (1964) Eq. 5.48, Primarily for debugging purposes
        aetup = [(self.semiMajorAxis[i], self.eccentricity[i]) for i in range(0, self.npoints)]
        def calc_c0(tuple_list):
            c0=[]
            for pair in tuple_list:
                a = pair[0]
                e = pair[1]
                c0i = a * (1-e**2) * (e**(-12/19)) * (1 + (121/304)*e**2 )**(-870/2299)
                c0.append(c0i)
            return c0
        c0 = calc_c0(aetup)
        if deviation_check:
            amplitude = max(c0) - min(c0)
            print(f"Calculated c0 values are within {amplitude:.4} of constant (maximum - minimum)")
        self.c0 = c0

    def set_points_per_orbit(self):
        df = np.diff(self.true_anomaly_deg)
        indices = np.where(np.sign(df) < 0)[0]
        extremes_t = [self.time[i] for i in indices]
        ppo = np.diff(indices)
        self.points_per_orbit = ppo
        self.points_per_orbit_avg = np.mean(ppo)

    def __init__(self, filename, i, j):
        
        self.data = load_spacehub_data(filename)
        self.i = i
        self.j = j
        print("load data complete")
        #set len and time
        print("Determining timesteps...")
        self.set_npoints()
        self.set_time()
        #set R vec
        print("Calculating orbital state vectors...")
        self.set_R_and_V()
        #Set mass attributes
        self.set_masses()
        #set hvec
        self.set_h_vector()
        #set N vec
        self.set_N_vector()
        #solve eccentricity vector
        self.set_ecc_vector()
        #solve Semi-major axis
        print("Calculating scalar orbital elements...")
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
        #Set orbit point-wise resolution
        self.set_points_per_orbit()
        #Set a-e relation constant c0 (Peters 1964 eq 5.48)
        self.set_c0()
        print("Done")

        

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

    def to_pandas(self):
        col_names = ["id", "x", "y", "z", "vx", "vy", "vz", "arg_of_peri_deg", "inclination_deg", "magN", "magR", "a", "true_anomaly_deg"]
        raise NotImplementedError("Not Implemented")
        pass

    def plot_orbit_3panel(self):
        df, i, j = self.data, self.i, self.j
        fig = plt.figure(figsize=(15, 5))
        ax0 = fig.add_subplot(131, projection='3d')
        ax1 = fig.add_subplot(132)
        ax2 = fig.add_subplot(133)
        #Plot axis spines
        ax0.quiver(0, 0, 0, 1, 0, 0, color='black', arrow_length_ratio=0.1, linewidth=2, label='X-axis')
        ax0.quiver(0, 0, 0,  0, 1, 0, color='black', arrow_length_ratio=0.1, linewidth=2, label='Y-axis')
        ax0.quiver(0, 0, 0, 0, 0, 1, color='black', arrow_length_ratio=0.1, linewidth=2, label='Z-axis')
        #plot orbit of j around i
        ax0.scatter(self.R_vec[0], self.R_vec[1], self.R_vec[2], s=0.3, c='red', zorder=2, label="Secondary")
        ax0.scatter(0, 0, 0, s=3, c='blue', zorder=1, label="Primary")
        #ax0.set_xlim(-8, 8)
        #ax0.set_ylim(-8, 8)
        ax1.scatter(0, 0, s=3, c='blue', zorder=1, label="Primary")
        ax1.scatter(self.R_vec[0], self.R_vec[1], s=0.3, c='red', zorder=2, label="Secondary")
        ax1.grid(visible=True, zorder=-1)
        ax1.set_xlim(-1,1)
        ax1.set_ylim(-1,1)

        ax2.scatter(self.R_vec[0], self.R_vec[2], s=0.3, c='red', zorder=2, label="Secondary")
        ax2.scatter(0, 0, s=3, c='blue', zorder=1, label="Primary")
        ax2.grid(visible=True, zorder=-1)

        return fig, [ax0, ax1, ax2]
    
    def plot_keplerian_evolution_basic(self, xlim = [], time_stop=0):
        #e, a, i, Omega, f, omega
        if len(xlim) < 2:
            xlim = [0, max(self.time)]

        fig, axes = plt.subplots(2,3, figsize=(15, 10))
        axs = axes.flatten()
        axs[0].plot(self.time, self.eccentricity)
        axs[0].set_ylim(-0.1, 1)
        #axs[0].set_ylabel("Eccentricity")
        axs[0].set_xlim(xlim)
        axs[0].set_title("Eccentricity")
        axs[0].set_ylabel(r"$e$")

        axs[1].plot(self.time, self.inclination_deg)
        axs[1].set_xlim(xlim)
        axs[1].set_ylim(0, 360)
        axs[1].set_title("Inclination")
        axs[1].set_ylabel(r"$i$ (deg)")

        axs[2].plot(self.time, self.magR)
        axs[2].set_xlim(xlim)
        axs[2].set_title("R (AU)")
        axs[2].set_ylim(0, 10)
        axs[2].set_title(r"Separation")
        axs[2].set_ylabel(r"$||R||$ ($AU$)")

        axs[3].plot(self.time, self.LongitudeAscendingNode_deg)
        axs[3].set_title("")
        axs[3].set_xlim(xlim)
        axs[3].set_ylim(0, 360)
        axs[3].set_ylabel(r"$\Omega$ (deg)")
        axs[3].set_title(r"Longitude of Ascending Node")

        axs[4].plot(self.time, self.semiMajorAxis)
        axs[4].set_xlim(xlim)
        #axs[4].set_ylim(0, 10)
        axs[4].set_ylabel(r"$a$ (AU)")
        axs[4].set_title(r"Semi-major Axis")
        axs[4].yaxis.set_major_formatter(FormatStrFormatter('%.2f'))
        
        axs[5].plot(self.time, self.argument_of_periapsis_deg)
        axs[5].set_xlim(xlim)
        axs[5].set_ylim(0, 360)
        axs[5].set_ylabel(r"$\omega$ (deg)")
        axs[5].set_title(r"Argument of Periapsis")

        for ax in axs:
            ax.set_xlabel(r"$yr (2\pi)^{-1}$")
        return fig, axs

    def plot_trajectory_3d(self, fig, ax, start_index = 0, *args, **kwargs):
        """
            args and kwargs are passed directly to animation.FuncAnimation
            pass fig, ax objects with 
        ``` fig = plt.figure()
            ax = fig.add_subplot(projection='3d')```
           Use HTML(ani.to_jshtml()) to render in IPython Notebooks
            
         """

        rx = self.R_vec[0][start_index:]
        ry = self.R_vec[1][start_index:]
        rz = self.R_vec[2][start_index:]
        global traj_points
        traj_points = ax.scatter3D(rx[0], ry[0], rz[0], c="purple")
        ref_frame_particle = ax.scatter3D(0, 0, 0, c="blue")
        traj = ax.plot(rx[0], ry[0], rz[0], c="gold")[0]


        def update(frame_num):
            # for each frame, update the data stored on each artist.
            window_size = 1
            window_start = frame_num-window_size if (frame_num > window_size) else 0
            x = rx[window_start:frame_num]
            y = ry[window_start:frame_num]
            z = rz[window_start:frame_num]
            # Update scatter
            data = np.stack([x,y,z]).T
            global traj_points  # Need global to reassign
            traj_points.remove()
            traj_points = ax.scatter3D(x, y, z, c="blue")
            # update the line plot:
            traj.set_data([rx[:frame_num], ry[:frame_num]])
            traj.set_3d_properties(rz[:frame_num])

            #traj.set_zdata(z[:frame])
            return (traj_points, traj)
        
        ani = animation.FuncAnimation(fig=fig, func=update, *args, **kwargs)
        return ani



class Theorize:



    def __init__(self):
        pass

    @staticmethod
    def enhancement_factor(e):
        sopra = 1 + (73/24)*e**2 + (37/96) * e**4
        sotto = (1-e**2)**(7/2)
        return sopra/sotto

    @staticmethod
    def decay_time_PN2p5(a0, m1, m2, e0):
        a0 = a0 * u.AU
        m1 = m1 * u.Msun
        m2 = m2 * u.Msun
        beta = (64/5) * G**3 * m1 * m2 * (m1+m2) * c**(-5) #Idk if the c^-5 is supposed to be here, but it fixes the units. G and c are not hard to write. stop using G=C=1
        f = Theorize.enhancement_factor(e0)

        T = a0**4 / (4*beta*f)
        return T.to(u.yr)

    @staticmethod
    def time_to_a_PN2p5(a_final):
        pass


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

def distance(data, key, i, j, npoints):
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
    for t in range(0, npoints):
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


#----Things I wrote and don't know what to do with but they might be useful at some point----
"""def a_theory(t): from peters 5.45 for a *circular orbit*
    t = t * u.yr / (np.pi * 2)
    m1 = orb.M_i * u.Msun
    m2 = orb.M_j * u.Msun
    beta = (64/5) * G**3 * m1 * m2 * (m1+m2) * c**(-5) #Idk if the c^-5 is supposed to be here, but it fixes the units. G and c are not hard to write. stop using G=C=1
    a0 = orb.semiMajorAxis[1] * u.AU
    a_theory = (a0**4 - 4*beta*t)**(1/4)
    return a_theory.value"""



#-------------Read/Write Operations--------------
# Load DefaultWriter output
def load_spacehub_data(filename, dropna=True):
    #units: AU = 1, year = 2pi, G = 1
    df = pd.read_csv(filename)
    if dropna:
        df = df.dropna()
    add_norms(df)
    
    return df

