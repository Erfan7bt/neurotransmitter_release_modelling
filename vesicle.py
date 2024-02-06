import numpy as np
import matplotlib.pyplot as plt

def poisson_spike_train(rate, t_sim, dt):
    survival = np.random.rand(int(t_sim/dt))
    spikes = np.where(survival < rate*dt, 1, 0)
    return spikes

def calcium_concentration(spike, tau_c, t_sim, dt):
    calcium = np.zeros(int(t_sim/dt))
    for i in range(1, len(spike)):
        calcium[i] = calcium[i-1]*(1-dt/tau_c) + spike[i-1]
    return calcium


def vesicle_fusion( tau_f, t_sim, dt, v_0, release_delay=0.1):
    release = np.zeros(int(t_sim/dt))
    release_delay = int(release_delay/dt)
    
    fusion = np.zeros(int(t_sim/dt))
    for i in range(1, len(release)):
        fusion[i] = dt/tau_f * (v_0 - release[i-1]) + (1 - dt/tau_f) * fusion[i-1]
        release[i] = calcium[i-release_delay]* spike[i-release_delay]* (v_0-fusion[i-release_delay])
    return fusion, release

rate = 10
t_sim = 10
dt = 0.001
tau_c = 0.4
tau_f = 0.3
v_0 = 10

spike = poisson_spike_train(rate, t_sim, dt)
calcium = calcium_concentration(spike, tau_c, t_sim, dt)
fusion, release = vesicle_fusion(tau_f, t_sim, dt, v_0)


plt.plot(spike, label='Spike', marker='o')
plt.plot(calcium, label='Calcium concentration')
plt.plot(fusion, label='Fusioned vesicle')
plt.plot(release, label='Release')
plt.xlabel('Time (ms)')
plt.ylabel('Spike')
plt.legend()
plt.title('Poisson spike train')
plt.savefig('vesicle.png')
plt.show()

