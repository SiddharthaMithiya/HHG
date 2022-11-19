""" Higher Harmonic Generation 1d using RK4 method """
import numpy as np
import matplotlib.pyplot as plt
import decorate_axes_dark as da
fig = plt.figure()
ax1 = fig.add_subplot(221)
ax3 = fig.add_subplot(223)
ax2 = fig.add_subplot(222)
ax4 = fig.add_subplot(224)
da.decorate_2d(fig, ax1)
da.decorate_2d(fig, ax2)
da.decorate_2d(fig, ax3)
da.decorate_2d(fig, ax4)

m_e = 9.1 * 10**-31         # kg
e = 1.6 * 10**-19           # C
c = 3 * 10**8               # m/s
ε0 = 8.85 * 10**-12         # F/m

λ = 800 * 10**-9           # m
I0 = 1 * 10**14 * 10**4     # Watt/m2
E0 = np.sqrt(2 * I0 / (c * ε0))
freq = c / λ                # Hz
Time_period = 1/freq        # sec

t0 = 0; tf = 10*Time_period; x0 = 0; xp0 = 0

ti = t0; xi = x0; xpi = xp0
t, dt = np.linspace(ti, tf, 500, retstep=True)
print(f't0 = {t0} sec ; t_f = {tf} sec\ninitial position = {x0} \ninitial velocity = {xp0} m/s\n')
print(f'I0 = {I0} \nE0 = {E0} \nλ = {round(λ * 10**9, 3)} nm')
print(f'Time period = {Time_period} sec \nfrequency = {freq} Hz')

def d2xdt2(t, x):
    # E_c = - (9 * 10**9) * e / x ** 2
    E_l = - (e / m_e) * E0 * np.sin(2*np.pi*freq*t) * np.sin(np.pi*freq*t/10) ** 2
    return E_l

E_l = E0 * np.sin(2*np.pi*freq*t)
E_pulse = -E0 * np.sin(2*np.pi*freq*t) * np.sin(np.pi*freq*t/10) ** 2
x = [xi]
for i in range(1, len(t)):
    J1 = dt * d2xdt2(ti, xi)
    K1 = dt * xpi
    J2 = dt * d2xdt2(ti + dt / 2, xi + K1 / 2)
    K2 = dt * (xpi + J1 / 2)
    J3 = dt * d2xdt2(ti + dt / 2, xi + K2 / 2)
    K3 = dt * (xpi + J2 / 2)
    J4 = dt * d2xdt2(ti + dt, xi + K3)
    K4 = dt * (xpi + J3)
    xpn = xpi + (J1 + 2 * J2 + 2 * J3 + J4) / 6
    xn = xi + (K1 + 2 * K2 + 2 * K3 + K4) / 6
    x.append(xn)
    ti += dt; xi = xn; xpi = xpn

x = np.array(x)                         # m
v = np.diff(x) / np.diff(t)             # m/s
Ek = (0.5 *m_e* v**2) * 10**19 / 1.6    # eV

#                           -:: Plotting ::-
fig.suptitle('Higher Harmonic Generation', color='yellow', fontsize=20)
ax1.plot(t, x)
ax2.plot(t[0: -1], Ek, c='m')
# ax3.plot(t, E_l,color='m')
ax3.plot(t, E_pulse, color='#83C167')
ax4.plot(Ek, x[0: -1])
ax3.axhline(0, color='white', lw=0.5)
ax1.set_xlabel(r'time $\longrightarrow$', fontsize=20)
ax1.set_ylabel(r'x $\longrightarrow$', fontsize=20)
ax2.set_xlabel(r'time $\longrightarrow$', fontsize=20)
ax2.set_ylabel(r'Ek $\longrightarrow$', fontsize=20)
ax3.set_xlabel(r'time $\longrightarrow$', fontsize=20)
ax3.set_ylabel(r'E laser $\longrightarrow$', fontsize=20)
ax4.set_xlabel(r'Ek $\longrightarrow$', fontsize=20)
ax4.set_ylabel(r'x $\longrightarrow$', fontsize=20)

fig.subplots_adjust(top=0.925, bottom=0.08, left=0.060, right=0.975, hspace=0.175, wspace=0.13)
plt.show()