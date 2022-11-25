""" Higher Harmonic Generation 2d using RK4 method """
import numpy as np
import matplotlib.pyplot as plt
plt.rc('font', **{'family': 'serif', 'size':15})

fig = plt.figure()
ax1 = fig.add_subplot(221)
ax4 = fig.add_subplot(222)
ax3 = fig.add_subplot(234)
ax2 = fig.add_subplot(235)
ax5 = fig.add_subplot(236)

# Defining the constants
m_e = 9.1 * 10**-31         # kg
e = 1.6 * 10**-19           # C
c = 3 * 10**8               # m/s
ε0 = 8.85 * 10**-12         # F/m

# Use λ = 800, 1200, 1600 nm
λ = 800 * 10**-9            # m
I0 = 1 * 10**14 * 10**4     # Watt/m2
E0 = np.sqrt(2 * I0 / (c * ε0))
freq = c / λ                # Hz
w = 2 * np.pi * freq        # Angular frequancy
Time_period = 1/freq        # sec
optical_cycle = 15          # (N) in the report
t0 = 0; tf = optical_cycle*Time_period

# Initial conditions
# x0 = 0; y0 = 0; xp0 = 0; yp0 = 0    # zero initial lateral velocity
x0 = 0; y0 = 0; xp0 = 60000; yp0 = 60000    # non zero initial lateral velocity


ti = t0; xi = x0; yi = y0; xpi = xp0; ypi = yp0
t, dt = np.linspace(ti, tf, 1000, retstep=True)
print(f't0 = {t0} sec ; t_f = {tf} sec\ninitial position = {x0} \ninitial velocity = {np.sqrt(xp0**2 + yp0**2)} m/s\n')
print(f'I0 = {I0} Watt/m2\nE0 = {E0} \nλ = {round(λ * 10**9, 3)} nm')
print(f'Time period = {Time_period} sec \nfrequency = {freq} Hz')
print('Ponderomotive energy = ', (e*E0/(2*m_e*w))**2 * 10 ** 19 / 1.6)
fig.suptitle('Higher Harmonic Generation', fontsize=20)

def d2xdt2(t, x):
    # Returns the acceleration of eletron
    acceleration = -(e / m_e) * E0 * np.sin(w * t) * np.sin(w * t / (2*optical_cycle)) ** 2
    return acceleration

def d2ydt2(t, y):
    # Returns the acceleration of eletron
    acceleration = -(e / m_e) * E0 * np.cos(w * t) * np.sin(w * t / (2*optical_cycle)) ** 2
    return acceleration


# This are the x and y component of electric field we used
E_pulse_x = E0 * np.sin(w * t) * np.sin(w * t / (2*optical_cycle)) ** 2
E_pulse_y = E0 * np.cos(w * t) * np.sin(w * t / (2*optical_cycle)) ** 2

# RK4 Method
x = [xi]; y = [yi]
for i in range(1, len(t)):
    J1_x = dt * d2xdt2(ti, xi)
    J1_y = dt * d2ydt2(ti, yi)

    K1_x = dt * xpi
    K1_y = dt * ypi

    J2_x = dt * d2xdt2(ti + dt / 2, xi + K1_x / 2)
    J2_y = dt * d2ydt2(ti + dt / 2, yi + K1_y / 2)

    K2_x = dt * (xpi + J1_x / 2)
    K2_y = dt * (ypi + J1_y / 2)

    J3_x = dt * d2xdt2(ti + dt / 2, xi + K2_x / 2)
    J3_y = dt * d2ydt2(ti + dt / 2, yi + K2_y / 2)

    K3_x = dt * (xpi + J2_x / 2)
    K3_y = dt * (ypi + J2_y / 2)

    J4_x = dt * d2xdt2(ti + dt, xi + K3_x)
    J4_y = dt * d2ydt2(ti + dt, yi + K3_y)

    K4_x = dt * (xpi + J3_x)
    K4_y = dt * (ypi + J3_y)

    xpn = xpi + (J1_x + 2 * J2_x + 2 * J3_x + J4_x) / 6
    ypn = ypi + (J1_y + 2 * J2_y + 2 * J3_y + J4_y) / 6

    xn = xi + (K1_x + 2 * K2_x + 2 * K3_x + K4_x) / 6
    yn = yi + (K1_y + 2 * K2_y + 2 * K3_y + K4_y) / 6

    x.append(xn)
    y.append(yn)

    ti += dt; xi = xn; xpi = xpn; yi = yn; ypi = ypn

x = np.array(x)                         # m
y = np.array(y)                         # m
r = np.sqrt(x**2 + y**2)                # m
v_x = np.diff(x) / np.diff(t)             # m/s
v_y = np.diff(y) / np.diff(t)             # m/s
Ek = (0.5 * m_e * (v_x ** 2 + v_y ** 2)) * 10 ** 19 / 1.6    # eV
print('Maximum kinetic energy : ', max(Ek), ' eV')





#                           -:: Plotting ::-
ax1.plot(x * 10**10, y * 10**10, color='red')
ax1.axis('equal')


ax2.plot(t[0: -1], Ek, c='green', lw=1.5)


ax3.plot(t, E_pulse_x, color='red', label='X component')
ax3.plot(t, E_pulse_y, color='blue', label='Y component')
ax3.legend(fontsize=13, framealpha=0.5, edgecolor='white')


ax4.plot(t, r * 10**10, color='blue')
ax4.axhline(4 * 0.529, linestyle='dashed', color='black')


ax5.plot(Ek, r[0: -1] * 10**10, color='red')
ax5.axhline(4 * 0.529, linestyle='dashed', color='black')


ax1.set_xlabel(r'x $\left( \AA \right) \longrightarrow$', fontsize=20)
ax1.set_ylabel(r'y $\left( \AA \right) \longrightarrow$', fontsize=20)
ax2.set_xlabel(r'time $\longrightarrow$', fontsize=20)
ax2.set_ylabel(r'Ek (eV) $\longrightarrow$', fontsize=20)
ax3.set_xlabel(r'time $\longrightarrow$', fontsize=20)
ax3.set_ylabel(r'E laser $\longrightarrow$', fontsize=20)
ax4.set_xlabel(r'time $\longrightarrow$', fontsize=20)
ax4.set_ylabel(r'r $\left( \AA \right) \longrightarrow$', fontsize=20)
ax5.set_xlabel(r'Ek (eV) $\longrightarrow$', fontsize=20)
ax5.set_ylabel(r'r $\left( \AA \right) \longrightarrow$', fontsize=20)



ax1.grid(True, lw=1, alpha=1, zorder=0)
ax2.grid(True, lw=1, alpha=1, zorder=0)
ax3.grid(True, lw=1, alpha=1, zorder=0)
ax4.grid(True, lw=1, alpha=1, zorder=0)
ax5.grid(True, lw=1, alpha=1, zorder=0)


fig.subplots_adjust(top=0.925, bottom=0.08, left=0.06, right=0.975, hspace=0.24, wspace=0.18)

plt.show()
