import numpy as np
from scipy.integrate import odeint
import matplotlib.pyplot as plt

# Hodgkin-Huxley parameters
gNa = 120.0  # Maximum conductances, mS/cm^2
gK = 36.0
gL = 0.3
ENa = 50.0   # Nernst reversal potentials, mV
EK = -77.0
EL = -54.387

# Time-dependent stimulus
def I_stim(t):
    return 10 * (10 <= t <= 20)

# Rate functions for the ion channels
def alpha_m(V): return 0.1*(V+40)/(1 - np.exp(-(V+40)/10))
def beta_m(V): return 4.0*np.exp(-(V+65)/18)
def alpha_h(V): return 0.07*np.exp(-(V+65)/20)
def beta_h(V): return 1/(1 + np.exp(-(V+35)/10))
def alpha_n(V): return 0.01*(V+55)/(1 - np.exp(-(V+55)/10))
def beta_n(V): return 0.125*np.exp(-(V+65)/80)

# Hodgkin-Huxley model equations
def hodgkin_huxley(X, t):
    V, m, h, n = X
    dVdt = (I_stim(t) - gNa*m**3*h*(V - ENa) - gK*n**4*(V - EK) - gL*(V - EL))
    dmdt = alpha_m(V)*(1-m) - beta_m(V)*m
    dhdt = alpha_h(V)*(1-h) - beta_h(V)*h
    dndt = alpha_n(V)*(1-n) - beta_n(V)*n
    return dVdt, dmdt, dhdt, dndt

# Initial conditions and time points
X0 = [-65, 0.05, 0.6, 0.32]
t = np.linspace(0, 50, 500)

# Solve ODE
X = odeint(hodgkin_huxley, X0, t)

# Plot results
plt.figure(figsize=(10, 8))
plt.subplot(2,1,1)
plt.title("Hodgkin-Huxley Neuron Model")
plt.plot(t, X[:,0], label="Membrane Potential (mV)")
plt.ylabel("Voltage (mV)")
plt.legend()

plt.subplot(2,1,2)
plt.plot(t, [I_stim(ti) for ti in t], label="Stimulus Current (uA/cm^2)")
plt.xlabel("Time (ms)")
plt.ylabel("Current (uA/cm^2)")
plt.legend()
plt.show()
