#!/usr/bin/env python3

import isa
from math import pi
import matplotlib.pyplot as plt

def sign(i):
	# My own sign function because no sign function exists in python
	# Ref: https://stackoverflow.com/questions/1986152/why-doesnt-python-have-a-sign-function
	if i<0:
		return -1
	else:
		return 1

# Simulating the ascent of a balloon

# Constants

g   = 9.8         # m/s²      - local gravitational acceleration
R   = 8.31        # J/(mol K) - gas constant
MHe = 0.004002602 # kg/mol    - helium molar mass

# Input data

m0 = 2.0 # kg - mass of balloon plus payload
D  = 2.0 # m³ - initial diameter of balloon
Cd = 0.6 #    - Drag coefficient
h0 = 0.0 # m  - Start height
Dp = 1.0 # m  - Diameter of parachute

# Calculations

p0  = isa.getPressure(h0)
T0  = isa.getTemperature(h0)
r0  = D/2              # m   - initial radius of balloon
V0  = 4/3 * pi * r0**3 # m³  - initial volume of balloon
n   = p0*V0/(R*T0)     # mol - amount of substance of Helium
mHe = n * MHe          # kg  - mass of helium

# Initial conditions

h  = 0   # height
v  = 0   # velocity
t  = 0   # time
dt = 0.1 # time step

m = m0 + mHe # kg - total mass
V = V0
h = h0

# Lists for plotting

tl   = []
hl   = []
vl   = []
Vl   = []
Dl   = []
rhol = []
pl   = []
Tl   = []

burst = False

while h>=0:

	T   = isa.getTemperature(h)
	p   = isa.getPressure(h)
	rho = isa.getDensity(h)

	tl.append(t/60)
	hl.append(h/1000)
	vl.append(v)
	pl.append(p)
	Tl.append(T)
	Dl.append(D)
	Vl.append(V)
	rhol.append(rho)

	t = t + dt

	if not burst and D>8:
		burst = True
		print("Burst time:     %.2f" % (t/60)   + " min")
		print("Burst altitude: %.2f" % (h/1000) + " km")

	if not burst:
		V = n*R*T/p # Volume of balloon is given by ideal gas law
		r = (3*V/(4*pi))**(1/3)
		D = 2*r
	else:
		V = 0 # Volume of balloon is zero
		D = 0
		r = Dp/2
		n = 0

	Fg = -m * g
	Fo = V * g * rho
	Fd = -0.5 * rho * v**2 * pi * r**2 * Cd * sign(v)

	F  = Fo + Fg + Fd

	a = F/m
	v = v + a*dt
	h = h + v*dt

print("Landing time:   %.2f" % (t/60)   + " min")

plt.subplot(231)
plt.plot(tl,hl)
plt.xlabel("Time / [min]")
plt.ylabel("Altitude / [km]")
plt.grid(True)

plt.subplot(232)
plt.plot(tl,vl)
plt.xlabel("Time / [min]")
plt.ylabel("Velocity / [m/s]")
plt.grid(True)

plt.subplot(233)
plt.plot(tl,Dl)
plt.xlabel("Time / [min]")
plt.ylabel("Diameter / [m]")
plt.grid(True)

plt.subplot(234)
plt.plot(pl,hl)
plt.xlabel("Pressure / [Pa]")
plt.ylabel("Altitude / [km]")
plt.grid(True)

plt.subplot(235)
plt.plot(Tl,hl)
plt.xlabel("Temperature / [K]")
plt.ylabel("Altitude / [km]")
plt.grid(True)

plt.subplot(236)
plt.plot(rhol,hl)
plt.xlabel("Density / [kg/m³]")
plt.ylabel("Altitude / [km]")
plt.grid(True)

plt.show()
