import numpy as np
import matplotlib.pyplot as plt

def poisson_spike_train(rate, t_sim, dt, silent=True): 
    survival = 1 - np.exp(-rate*dt) 
    spike = np.random.rand(int(t_sim/dt)) < survival
    if silent:
        #kill the last 10 percent of the spikes
        spike[-int(len(spike)*0.1):] = 0
    return spike

def calcium_concentration(spike, tau_c, t_sim, dt, R=10):
    spike = spike * R #R is the amount of calcium ions released per spike
    calcium = np.zeros(int(t_sim/dt))
    for i in range(1, len(spike)):
        calcium[i] = calcium[i-1]*(1-dt/tau_c) + dt/tau_c* spike[i-1]
    return calcium

def vesicle_fusion( tau_f, t_sim, dt, v_0, release_delay):
    release = np.zeros(int(t_sim/dt))
    release_delay = int(release_delay/dt)
    
    fusion = np.zeros(int(t_sim/dt))
    for i in range(1, len(release)):
        fusion[i] = dt/tau_f * (v_0 - release[i-1]) + (1 - dt/tau_f) * fusion[i-1]
        release[i] = calcium[i-release_delay]* spike[i-release_delay]* (v_0-fusion[i-release_delay])
    return fusion, release

rate = 30
t_sim = 1
dt = 0.001
tau_c = 30 *dt
tau_f = 100 *dt #time constant for vesicle fusion process higher than calcium influx time constant
v_0 = 30 #maximum vesicle pool size
release_delay = 0.03 #delay between calcium influx and vesicle release

spike = poisson_spike_train(rate, t_sim, dt)
calcium = calcium_concentration(spike, tau_c, t_sim, dt)
fusion, release = vesicle_fusion(tau_f, t_sim, dt, v_0, release_delay)


fig, ax = plt.subplots(3,1, figsize=(10,10), layout='constrained')
ax[0].plot(spike)
ax[0].set_title(f'Spike Train with {rate} Hz rate')
ax[0].set_xlabel('Time (s)')

ax[1].plot(calcium)
ax[1].set_title('Calcium Concentration in the presynaptic terminal')
ax[1].set_ylabel('Concentration')
ax[1].set_xlabel('Time (s)')

ax[2].plot(fusion, label='Fused Vesicles')
ax[2].plot(release, label='Released neurotransmitter')
ax[2].set_xlabel('Time (s)')
ax[2].set_ylabel('#')
ax[2].set_title('Vesicle Pool')
ax[2].legend()
plt.show()

class Vesicle:

    def __init__(self, rate, t_sim, dt, tau_c, tau_f, v_0, release_delay):
        self.rate = rate
        self.t_sim = t_sim
        self.dt = dt
        self.tau_c = tau_c
        self.tau_f = tau_f
        self.v_0 = v_0
        self.release_delay = release_delay
        self.spike = poisson_spike_train(self.rate, self.t_sim, self.dt)
        self.calcium = calcium_concentration(self.spike, self.tau_c, self.t_sim, self.dt)
        self.fusion, self.release = vesicle_fusion(self.tau_f, self.t_sim, self.dt, self.v_0, self.release_delay)

#trials 
trials = 10
for i in range(trials):
    v = Vesicle(rate, t_sim, dt, tau_c, tau_f, v_0, release_delay)

avg_spike = np.mean(v.spike) 
avg_fusion = np.mean(v.fusion)
avg_release = np.mean(v.release)

