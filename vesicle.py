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


def vesicle_fusion( tau_f, t_sim, dt, v_0):
    release = np.zeros(int(t_sim/dt))
    
    fusion = np.zeros(int(t_sim/dt))
    for i in range(1, len(release)):
        fusion[i] = dt/tau_f * (v_0 - release[i-1]) + (1 - dt/tau_f) * fusion[i-1]
        release[i] = calcium[i]* spike[i]* fusion[i]
    return fusion


spike = poisson_spike_train(17, 1, 0.001)

calcium = calcium_concentration(spike, 0.4, 1, 0.001)

fusion = vesicle_fusion(0.1, 1, 0.001, 10)
print(fusion)

plt.plot(spike, label='Spike')
plt.plot(calcium, label='Calcium')
plt.plot(fusion, label='Fusion')
plt.xlabel('Time (ms)')
plt.ylabel('Spike')
plt.legend()
plt.title('Poisson spike train')
plt.show()

