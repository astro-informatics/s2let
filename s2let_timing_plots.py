
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D

fontsize = 13
markers = []
for m in Line2D.markers:
    try:
        if len(m) == 1 and m != ' ':
            markers.append(m)
    except TypeError:
        pass
styles = markers + [
    r'$\lambda$',
    r'$\bowtie$',
    r'$\circlearrowleft$',
    r'$\clubsuit$',
    r'$\checkmark$']

fig, axs = plt.subplots(1,2, figsize=(12,4.5))
axs = axs.ravel()

suffix = 'L128_s0_B3_N5'
data_L_N = np.genfromtxt('timings_errors_'+suffix+'.csv', delimiter=';', names=True)
N = int(data_L_N['N'][0])
J = int(data_L_N['J'][0])
J_min = int(data_L_N['J_min'][0])
spin = int(data_L_N['spin'][0])
B = double(data_L_N['B'][0]) #DOUBLE
outname = 's='+str(spin)+', N='+str(N)+', B='+str(B)+', Jmin='+str(J_min)+', J='+str(J) 

ind1 = (data_L_N['multires'] == 1)
ind2 = (data_L_N['multires'] == 0)
L = data_L_N['L'][ind1]
p0 = axs[0].plot(L, 2e-7 * L**3, color='red')
p1 = axs[0].plot(data_L_N['L'][ind1], data_L_N['min_duration_inverse'][ind1], color='blue', ls='solid', marker=styles[0])
p2 = axs[0].plot(data_L_N['L'][ind1], data_L_N['min_duration_forward'][ind1], color='black', ls='solid', marker=styles[1])
p3 = axs[0].plot(data_L_N['L'][ind2], data_L_N['min_duration_inverse'][ind2], color='blue', ls='dashed', marker=styles[0])
p4 = axs[0].plot(data_L_N['L'][ind2], data_L_N['min_duration_forward'][ind2], color='black', ls='dashed', marker=styles[1])
lg = axs[0].legend([p1[0],p2[0],p3[0],p4[0],p0[0]],
              ['multi resolution, min_duration_inverse', 'multi resolution, min_duration_forward',
                'full resolution, min_duration_inverse', 'full resolution, min_duration_forward', 'L^3 scaling'],
                loc='upper left', fontsize=fontsize)
lg.draw_frame(False)
axs[0].set_yscale('log')
axs[0].set_xscale('log')
axs[0].set_ylim([1e-4, 1e4])
axs[0].set_xlim([data_L_N['L'].min()/1.5, data_L_N['L'].max()*1.5])
axs[0].set_xlabel('L', fontsize=fontsize)
axs[0].set_ylabel('Average duration', fontsize=fontsize)
axs[0].set_title('S2LET timing scaling: '+outname, fontsize=fontsize)

p7 = axs[1].plot(L, 1e-15 * L, color='red')
p5 = axs[1].plot(data_L_N['L'][ind1], data_L_N['avg_error'][ind1], color='blue', ls='solid', marker=styles[0])
p6 = axs[1].plot(data_L_N['L'][ind2], data_L_N['avg_error'][ind2], color='black', ls='dashed', marker=styles[0])
lg = axs[1].legend([p5[0],p6[0],p7[0]],
              ['multi resolution, avg_error', 'full resolution, avg_error', 'L scaling'],
                loc='upper left', fontsize=fontsize)
lg.draw_frame(False)
axs[1].set_yscale('log')
axs[1].set_xscale('log')
axs[1].set_xlim([data_L_N['L'].min()/1.5, data_L_N['L'].max()*1.5])
axs[1].set_xlabel('L', fontsize=fontsize)
axs[1].set_ylabel('Average error', fontsize=fontsize)
axs[1].set_title('S2LET error scaling: '+outname, fontsize=fontsize)
axs[1].set_ylim([1e-15, 1e-12])


fig.tight_layout()
fig.savefig('s2let_timing_'+suffix+'.png')
plt.show()



