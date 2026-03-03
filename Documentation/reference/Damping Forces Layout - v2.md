**→Both of the approaches below derive the radial component. Theoretically if everything is self consistent, the $\bar{R}$ derived from both will be the same (I think)**
### Derivation from semi-major axis evolution

Use MD 2.149 in combo with Fairbairn definitions
$$
\begin{array}{l}
\frac{1}{m}\frac{dL}{dt}=r\bar{T}; \\ 
\tau_L^{-1}=\frac{T_{net}}{L};\\
\frac{dL}{dt}=-T_{net};\\
\dot{a}=\frac{da}{dt}=-a\tau_a^{-1}
\end{array}
$$
From the datacube, $T_{net}$ is known, and thus $\bar{T}$ can be calculated with a known $R$ and $L$ 
$$\bar{T}=-\frac{T_{net}}{mr}$$
Substitute this and the $\dot{a}$ relation into equation 2.145
$$\frac{da}{dt}=2\frac{a^{3/2}}{\sqrt{\mu(1-e^2)}}[\bar{R}e\sin{f} + \bar{T}(1+e\cos{f})]=2\frac{a^{3/2}}{\sqrt{\mu(1-e^2)}}[\bar{R}e\sin{f} + -\frac{T_{net}}{mr}(1+e\cos{f})]$$
$$\dot{a}=-a\tau_a^{-1}=2\frac{a^{3/2}}{\sqrt{\mu(1-e^2)}}[\bar{R}e\sin{f} + -\frac{T_{net}}{mr}(1+e\cos{f})]$$
Solve for $\bar{R}$
$$\boxed{\bar{R}=\frac{1}{e\sin f}\left[\frac{-a\tau_a^{-1}\sqrt{\mu(1-e^2)}}{2a^{3/2}}+\frac{T_{net}}{mr}(1+e\cos f)\right]}$$
**NOTE: This is undefined for circular orbits, but we still expect migration torques in the circular case...?**

Substitutions for implementation into SpaceHub
$$h=|\vec{h}|=|\vec{r}\times \vec{v}|=\sqrt{\mu a (1-e^2)}$$
$$e\sin f = \frac{h\dot{r}}{\mu}$$
$$1+e\cos f = \frac{h^2}{\mu r}$$
$$\cos E = \frac{1}{e}\left(1-\frac{r}{a}\right) \quad (\mathrm{2.42})$$
Final acceleration vector for $e, a$ damping forces
$$\vec{a}_{migration}=\bar{R}\hat{r}+\bar{T}\hat{\theta}+0\hat{z}=\bar{R}\hat{r}+\bar{T}\left(\frac{\hat{z}\times\vec{r}}{r}\right)+0\hat{z}$$

### Derivation from eccentricity damping 
Use:
$$
\begin{array}{l}
\frac{dh}{dt}=\frac{1}{m}\frac{dL}{dt}=r\bar{T}; \\ 
\tau_L^{-1}=\frac{T_{net}}{L};\\
\frac{dL}{dt}=-T_{net};\\
\dot{e}=\frac{de}{dt}=-e\tau_e^{-1}
\end{array}
$$
$$\frac{de}{dt}=\sqrt{a\mu^{-1}(1-e^2)}\left[\bar{R}\sin f + \bar{T} (\cos f + \cos E) \right]$$
Where $E$ here is the eccentric anomaly, not the energy. Substitute in known $\bar{T}$ and $\dot{e}$
$$\dot{e}=-e\tau_e^{-1}=\sqrt{a\mu^{-1}(1-e^2)}\left[\bar{R}\sin f - \frac{T_{net}}{mr} (\cos f + \cos E) \right]$$
Solve for radial component $\bar{R}$

$$\frac{-e\tau_e^{-1}}{\sqrt{a\mu^{-1}(1-e^2)}}+\frac{T_{net}}{mr}(\cos f + \cos E)=\bar{R}\sin f$$

$$\bar{R}=\frac{1}{\sin f}\left(\frac{-e\tau_e^{-1}}{\sqrt{a\mu^{-1}(1-e^2)}}+\frac{T_{net}}{mr}(\cos f + \cos E)\right)$$



### SpaceHub-friendly process

Particle position $\vec{p}$  and velocity $\vec{v}$ relative to the COM frame which the disk is at rest in
Disk parameters evaluated at particle position
##### Calculate $p$ and $q$ to access datacube
Fairbairn disk prescription:
$$\Sigma(r)=\Sigma_p(r/a)^{-p}$$
$$c_s(r)=c_{s,p}(r/a)^{-q/2}$$
$$P=c_s^2\Sigma\quad \mathrm{(EOS)}$$
$$T_{disc}(r)\propto\frac{P}{\Sigma}=c_s^2\propto r^{-q}$$
$$\mathrm{aspect\ ratio\ } h=\frac{H(r)}{r}=h_p(r/a)^{(1-q)/2}$$
We do not care about $\Sigma_p, h_p,$ etc because $p$ and $q$ depend only on the local log gradient (?)

from EOS
$$\ln P = \ln(c_s^2)+\ln(\Sigma)$$
$$\frac{\ln P}{d\ln r} = \frac{\ln(c_s^2)}{d\ln r}+\frac{\ln(\Sigma)}{d\ln r}$$
Substitute in derivatives based on power law definitions
$$\frac{d\ln P}{d\ln r}=-q-p$$
$$p=-\frac{d\ln\Sigma}{d\ln r}, \quad q = \frac{dln\Sigma}{d\ln r}-\frac{dlnP}{d\ln r}$$
Then extract from datacube (where $h$ is the aspect ratio $h=H(r)/r$
$$\tau_a^{-1}(p,q,h,e), \quad T_{net}(p,q,h,e)$$
#### Calculate orbital elements from position/velocity

$$h=|\vec{h}|=|\vec{r}\times \vec{v}|=\sqrt{\mu a (1-e^2)}$$
radial velocity from dot product of position and velocity vectors:
$$\dot{r}=\frac{\vec{r}\cdot\vec{v}}{|\vec{r}|}$$
$$e\sin f = \frac{h\dot{r}}{\mu}$$
$$1+e\cos f = \frac{h^2}{\mu r}$$
$$\cos E = \frac{1}{e}\left(1-\frac{r}{a}\right) \quad (\mathrm{2.42})$$
Calculate force vector components with the above:

$$\bar{T}=-\frac{T_{net}}{mr}$$
$$\bar{R}=\frac{1}{e\sin f}\left[\frac{-a\tau_a^{-1}\sqrt{\mu(1-e^2)}}{2a^{3/2}}+\frac{T_{net}}{mr}(1+e\cos f)\right]$$
$$\bar{N}=0$$
Then convert $\langle\bar{R}, \bar{T}, 0\rangle$ to cartesian and add to the acceleration of the particle (and subtract from the central mass, if using assumption that disk is anchored to central mass and not at rest in COM frame)