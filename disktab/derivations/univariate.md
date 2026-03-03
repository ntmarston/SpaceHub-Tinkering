# Univariate Derivations

## Definitions

$$f=1-\sqrt{\frac{R_*}{R}}$$

$$\Omega=\sqrt{\frac{GM}{R^3}}$$

## Equation I

$$\frac{k_B T_c}{\mu m_p} + \frac{4\sigma T_c^4}{3c\rho}-\left[\frac{\dot{M}\Omega^2f}{3\pi\alpha\rho}\right]^{2/3}=0$$

$$A+B\rho^{-1}=C\rho^{-2/3}$$

where:

$$A=\frac{k_B T_c}{\mu m_p}$$

$$B=\frac{4\sigma T_c^4}{3c}$$

$$C=\left[\frac{f\dot{M}\Omega^2}{3\pi\alpha}\right]^{2/3}$$

## Equation II

$$\rho^{2/3}(\kappa_0+\kappa_1\rho T_c^{-7/2})-\left[\frac{32\pi\sigma\Omega R^3}{9GM\dot{M}f}\right]T_c^4\left[\frac{3\pi\alpha}{\dot{M}f\Omega^2}\right]^{1/3}=0$$

$$\rho^{2/3}(D+E\rho)=F$$

where:

$$D=\kappa_0$$

$$E=\kappa_1 T_c^{-7/2}$$

$$F=\left[\frac{32\pi\sigma\Omega R^3}{9GM\dot{M}f}\right]T_c^4\left[\frac{3\pi\alpha}{\dot{M}f\Omega^2}\right]^{1/3}$$

## Quadratic Solution

$$\left(\frac{CE}{F}\right)\rho^2+\left(\frac{CD}{F}-A\right)\rho-B=0$$

with coefficients:
- $a=\left(\frac{CE}{F}\right)$
- $b=\left(\frac{CD}{F}-A\right)$
- $c=-B$

$$\rho=\frac{-\left(\frac{CD}{F}-A\right)\pm\sqrt{\left(\frac{CD}{F}-A\right)^2-4\left(\frac{CE}{F}\right)(-B)}}{2\left(\frac{CE}{F}\right)}$$

## Temperature Derivatives

Writing coefficients separately from their $T_c$ dependence:

**Function $f(T_c)$:**

$$\rho^{2/3}\left(\frac{E\rho}{T_c^{7/2}}+D\right)-FT_c^4$$

**Derivative $f'(T_c)$:**

$$\frac{E \rho\left(T_{c}\right) \left(10T_{c} \frac{d\rho}{dT_c}\left(T_{c}\right) - 21 \rho\left(T_{c}\right)\right) + 4T_{c}^{\frac{9}{2}} \left(D \frac{d\rho}{dT_c}\left(T_{c}\right) - 6FT_{c}^{3} \sqrt[3]{\rho\left(T_{c}\right)}\right)}{6T_{c}^{\frac{9}{2}} \sqrt[3]{\rho\left(T_{c}\right)}}$$
