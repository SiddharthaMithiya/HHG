""" RK4 ; 2nd order differential equation """
import numpy as np
import matplotlib.pyplot as plt
plt.rc('font', **{'family': 'serif', 'size':15})

fig = plt.figure()
ax1 = fig.add_subplot(221)
ax3 = fig.add_subplot(223)
ax2 = fig.add_subplot(222)
ax4 = fig.add_subplot(224)

t0 = 0; tf = 1; x0 = 00; xp0 = 0
ti = t0; xi = x0; xpi = xp0
t, dt = np.linspace(ti, tf, 500, retstep=True)
print(f't0 = {t0} ; t_f = {tf} \ninitial position = {x0} \ninitial velocity={xp0}')

def d2xdt2(t, x):
    E_l = 0.5*np.sin(2*np.pi*2*t)
    return E_l

E_l = 0.5*np.sin(2*np.pi*2*t)

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

x = np.array(x)
v = np.diff(x) / np.diff(t)
Ek = 0.5 * v**2

#                           -:: Plotting ::-

ax1.plot(t, x)
ax2.plot(t[0: -1], Ek, c='m')
ax3.plot(t, E_l, color='red')
ax4.plot(x[0: -1], Ek)


ax1.grid(True, lw=1, alpha=1, zorder=0)
ax2.grid(True, lw=1, alpha=1, zorder=0)
ax3.grid(True, lw=1, alpha=1, zorder=0)
ax4.grid(True, lw=1, alpha=1, zorder=0)


ax3.axhline(0, color='white', lw=0.5)
ax1.set_xlabel(r'time $\longrightarrow$', fontsize=20)
ax1.set_ylabel(r'x $\longrightarrow$', fontsize=20)
ax2.set_xlabel(r'time $\longrightarrow$', fontsize=20)
ax2.set_ylabel(r'Ek $\longrightarrow$', fontsize=20)
ax3.set_xlabel(r'time $\longrightarrow$', fontsize=20)
ax3.set_ylabel(r'E laser $\longrightarrow$', fontsize=20)
ax4.set_xlabel(r'x $\longrightarrow$', fontsize=20)
ax4.set_ylabel(r'Ek $\longrightarrow$', fontsize=20)

fig.subplots_adjust(top=0.925, bottom=0.08, left=0.060, right=0.975, hspace=0.175, wspace=0.185)
plt.show()
