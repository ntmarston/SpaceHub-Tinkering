# Analytical Logarithmic Gradients for the Sirko–Goodman Self-Regulating Disk

## Setup

In the self-regulating (Sirko & Goodman) zone, the Toomre parameter is fixed at $Q = Q_{\min}$. The governing equations are:

$$\rho = \frac{\Omega^2}{2\pi G Q_{\min}} \tag{1}$$

$$c_s = \left(\frac{\dot{M} f \Omega^2}{6\pi\alpha\rho}\right)^{1/3} \tag{2}$$

$$P = c_s^2\,\rho \tag{3}$$

$$H = \frac{c_s}{\Omega} \tag{4}$$

$$\nu = \alpha\, c_s\, H \tag{5}$$

$$\Sigma = \frac{\dot{M} f}{3\pi\nu} = \frac{\dot{M} f\, \Omega}{3\pi\alpha\, c_s^2} \tag{6}$$

where $\Omega = \sqrt{GM/R^3}$ and $f = 1 - \sqrt{R_*/R}$.

The central temperature $T_c$ is determined by the energy balance:

$$\mathcal{C}_1\, T_c^4 + \mathcal{C}_2\, T_c = \mathcal{C}_3 \tag{7}$$

with coefficients:

$$\mathcal{C}_1 = \frac{4\sigma}{3c} \quad\text{(constant — no $R$ dependence)} \tag{8}$$

$$\mathcal{C}_2 = \frac{\rho\, k_B}{\mu\, m_p} \tag{9}$$

$$\mathcal{C}_3 = \rho\left(\frac{\dot{M} f \Omega^2}{6\pi\alpha\rho}\right)^{2/3} = \rho^{1/3}\left(\frac{\dot{M} f \Omega^2}{6\pi\alpha}\right)^{2/3} \tag{10}$$

**Key observation:** In this regime, $\rho$, $c_s$, $P$, $H$, $\nu$, and $\Sigma$ are all **explicit** algebraic functions of $R$ (and global constants $M, \dot{M}, \alpha$, etc.). None of them depend on $T_c$. Only $T_c$ itself requires solving the implicit equation (7).

We define the stored gradients as:

$$\texttt{grad\_X} \equiv -\frac{d\ln X}{d\ln R}$$

## Preliminary: Logarithmic Derivatives of Auxiliary Quantities

$$\frac{d\ln\Omega}{d\ln R} = -\frac{3}{2} \tag{11}$$

$$\frac{d\ln\rho}{d\ln R} = 2\,\frac{d\ln\Omega}{d\ln R} = -3 \tag{12}$$

For $f = 1 - \xi$ where $\xi \equiv \sqrt{R_*/R}$:

$$\frac{df}{dR} = \frac{1}{2}\sqrt{\frac{R_*}{R^3}} = \frac{\xi}{2R}$$

$$\frac{d\ln f}{d\ln R} = \frac{R}{f}\cdot\frac{df}{dR} = \frac{\xi}{2f} \tag{13}$$

## Gradient of $c_s$

From equation (2), $c_s^3 = \dot{M} f \Omega^2 / (6\pi\alpha\rho)$, so:

$$3\,\frac{d\ln c_s}{d\ln R} = \frac{d\ln f}{d\ln R} + 2\,\frac{d\ln\Omega}{d\ln R} - \frac{d\ln\rho}{d\ln R}$$

$$= \frac{\xi}{2f} + 2\!\left(-\frac{3}{2}\right) - (-3)$$

$$= \frac{\xi}{2f} - 3 + 3 = \frac{\xi}{2f}$$

$$\boxed{\frac{d\ln c_s}{d\ln R} = \frac{\xi}{6f}} \tag{14}$$

## Gradient of $P$ (Pressure)

From $P = c_s^2\,\rho$:

$$\frac{d\ln P}{d\ln R} = 2\,\frac{d\ln c_s}{d\ln R} + \frac{d\ln\rho}{d\ln R}$$

$$= \frac{2\xi}{6f} + (-3) = \frac{\xi}{3f} - 3$$

$$\boxed{\texttt{grad\_P} = -\frac{d\ln P}{d\ln R} = 3 - \frac{\xi}{3f}} \tag{15}$$

**Check:** For $R \gg R_*$, $\xi \to 0$, $f \to 1$, so $\texttt{grad\_P} \to 3$.

## Gradient of $\Sigma$ (Surface Density)

From equation (6), $\Sigma = \dot{M} f\, \Omega / (3\pi\alpha\, c_s^2)$:

$$\frac{d\ln\Sigma}{d\ln R} = \frac{d\ln f}{d\ln R} + \frac{d\ln\Omega}{d\ln R} - 2\,\frac{d\ln c_s}{d\ln R}$$

$$= \frac{\xi}{2f} + \left(-\frac{3}{2}\right) - \frac{2\xi}{6f}$$

$$= \frac{\xi}{2f} - \frac{\xi}{3f} - \frac{3}{2}$$

$$= \frac{3\xi - 2\xi}{6f} - \frac{3}{2} = \frac{\xi}{6f} - \frac{3}{2}$$

$$\boxed{\texttt{grad\_Sigma} = -\frac{d\ln\Sigma}{d\ln R} = \frac{3}{2} - \frac{\xi}{6f}} \tag{16}$$

**Check:** For $R \gg R_*$, $\texttt{grad\_Sigma} \to 3/2$.

## Gradient of $T_c$ (Central Temperature)

This is the only quantity that requires implicit differentiation. Start from equation (7):

$$\mathcal{C}_1\, T_c^4 + \mathcal{C}_2\, T_c = \mathcal{C}_3$$

First, compute the logarithmic derivatives of the coefficients.

### $\mathcal{C}_1$

Constant $\Rightarrow$ $\dfrac{d\ln\mathcal{C}_1}{d\ln R} = 0$.

### $\mathcal{C}_2$

From equation (9), $\mathcal{C}_2 = \rho\, k_B/(\mu m_p) \propto \rho$, so:

$$\frac{d\ln\mathcal{C}_2}{d\ln R} = \frac{d\ln\rho}{d\ln R} = -3 \tag{17}$$

### $\mathcal{C}_3$

From equation (10), $\mathcal{C}_3 = \rho^{1/3}\bigl(\dot{M} f \Omega^2/(6\pi\alpha)\bigr)^{2/3}$, so:

$$\frac{d\ln\mathcal{C}_3}{d\ln R} = \frac{1}{3}\,\frac{d\ln\rho}{d\ln R} + \frac{2}{3}\!\left[\frac{d\ln f}{d\ln R} + 2\,\frac{d\ln\Omega}{d\ln R}\right]$$

$$= \frac{1}{3}(-3) + \frac{2}{3}\!\left[\frac{\xi}{2f} + 2\!\left(-\frac{3}{2}\right)\right]$$

$$= -1 + \frac{2}{3}\!\left[\frac{\xi}{2f} - 3\right]$$

$$= -1 + \frac{\xi}{3f} - 2 = -3 + \frac{\xi}{3f} \tag{18}$$

### Implicit differentiation

Taking $d/d(\ln R)$ of both sides of equation (7):

$$4\,\mathcal{C}_1\, T_c^4 \cdot g + \mathcal{C}_2\, T_c\cdot(-3) + \mathcal{C}_2\, T_c \cdot g = \mathcal{C}_3\!\left(-3 + \frac{\xi}{3f}\right) \tag{19)$$

where $g \equiv \dfrac{d\ln T_c}{d\ln R}$.

Solving for $g$:

$$g\,\bigl(4\,\mathcal{C}_1\, T_c^4 + \mathcal{C}_2\, T_c\bigr) = \mathcal{C}_3\!\left(-3 + \frac{\xi}{3f}\right) + 3\,\mathcal{C}_2\, T_c$$

$$\boxed{\texttt{grad\_T} = -g = -\frac{\mathcal{C}_3\!\left(-3 + \dfrac{\xi}{3f}\right) + 3\,\mathcal{C}_2\, T_c}{4\,\mathcal{C}_1\, T_c^4 + \mathcal{C}_2\, T_c}} \tag{20}$$

### Limiting cases

**Radiation-dominated** ($\mathcal{C}_1 T_c^4 \gg \mathcal{C}_2 T_c$, i.e. $\mathcal{C}_3 \approx \mathcal{C}_1 T_c^4$), with $\xi \to 0$:

$$g \approx \frac{\mathcal{C}_1 T_c^4\cdot(-3) + 0}{4\,\mathcal{C}_1 T_c^4} = -\frac{3}{4} \quad\Rightarrow\quad \texttt{grad\_T} \to \frac{3}{4}$$

**Gas-dominated** ($\mathcal{C}_2 T_c \gg \mathcal{C}_1 T_c^4$, i.e. $\mathcal{C}_3 \approx \mathcal{C}_2 T_c$), with $\xi \to 0$:

$$g \approx \frac{\mathcal{C}_2 T_c\cdot(-3) + 3\,\mathcal{C}_2 T_c}{\mathcal{C}_2 T_c} = 0 \quad\Rightarrow\quad \texttt{grad\_T} \to 0$$

## Summary

All gradients for the self-regulating zone, with $\xi = \sqrt{R_*/R}$ and $f = 1 - \xi$:

| **Gradient** | **Formula** | **Limit ($R \gg R_*$)** |
|---|---|---|
| $\texttt{grad\_T} = -\dfrac{d\ln T_c}{d\ln R}$ | Eq. (20) (implicit) | $3/4$ |
| $\texttt{grad\_Sigma} = -\dfrac{d\ln\Sigma}{d\ln R}$ | $\dfrac{3}{2} - \dfrac{\xi}{6f}$ | $3/2$ |
| $\texttt{grad\_P} = -\dfrac{d\ln P}{d\ln R}$ | $3 - \dfrac{\xi}{3f}$ | $3$ |

**Note:** $\texttt{grad\_P} = 2\,\texttt{grad\_Sigma}$ for all $R$, which follows from $P = c_s^2\rho$ and $\Sigma \propto f\Omega/c_s^2$, giving $P/\Sigma \propto c_s^4 \rho / (f\Omega) \propto \Sigma$ after simplification. This provides a useful cross-check.
