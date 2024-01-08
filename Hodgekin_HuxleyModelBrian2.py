from brian2 import *

# Hodgkin-Huxley模型参数
area = 20000*umetre**2
Cm = 1*ufarad*cm**-2 * area
gNa = 120*msiemens*cm**-2 * area
gK = 36*msiemens*cm**-2 * area
gL = 0.3*msiemens*cm**-2 * area
ENa = 50*mV
EK = -77*mV
EL = -54.4*mV

# Hodgkin-Huxley模型方程
eqs = '''
dv/dt = (gNa*m**3*h*(ENa-v) + gK*n**4*(EK-v) + gL*(EL-v) + I)/Cm : volt
dm/dt = alpha_m*(1-m) - beta_m*m : 1
dn/dt = alpha_n*(1-n) - beta_n*n : 1
dh/dt = alpha_h*(1-h) - beta_h*h : 1
alpha_m = 0.1*(mV**-1)*(25*mV-v+VT)/
  (exp((25*mV-v+VT)/(10*mV)) - 1)/ms : Hz
beta_m = 4*exp((v-VT-60*mV)/(18*mV))/ms : Hz
alpha_h = 0.07*exp((v-VT-50*mV)/(20*mV))/ms : Hz
beta_h = 1/(exp((30*mV-v+VT)/(10*mV)) + 1)/ms : Hz
alpha_n = 0.01*(mV**-1)*(10*mV-v+VT)/
  (exp((10*mV-v+VT)/(10*mV)) - 1)/ms : Hz
beta_n = 0.125*exp((v-VT-60*mV)/(80*mV))/ms : Hz
I : amp
'''

# 创建一个神经元，运行模拟
neuron = NeuronGroup(1, eqs, method='exponential_euler')
neuron.v = -65*mV
neuron.I = '1*nA'
run(20*ms)

# 绘图
plot(neuron.v)
xlabel('Time (ms)')
ylabel('Membrane potential (mV)')
show()
