from __future__ import print_function
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

Lminplot = 20
Lplot = 2**np.arange(5, 12)

for suffix in ['L2048_s0_B2_N5', 'L2048_s2_B2_N5']:

    fig, axs = plt.subplots(1, 1, figsize=(6, 3.0))
    # axs = axs.ravel()
    data_L_N = np.genfromtxt('timings_errors_' + suffix + '.csv', delimiter=';', names=True)
    N = int(data_L_N['N'][0])
    J = int(data_L_N['J'][0])
    J_min = int(data_L_N['J_min'][0])
    spin = int(data_L_N['spin'][0])
    B = int(data_L_N['B'][0])
    outname = 's=' + str(spin) + ', N=' + str(N) + ', B=' + str(B) + ', Jmin=' + str(J_min)

    ind1 = np.logical_and(data_L_N['multires'] == 1, data_L_N['L'] > Lminplot)
    ind2 = np.logical_and(data_L_N['multires'] == 0, data_L_N['L'] > Lminplot)
    L = data_L_N['L'][ind1]
    Lbis = data_L_N['L'][ind1]
    print(Lbis)
    p0 = axs.plot(Lbis, 1e-8 * Lbis**3, color='red')
    p1 = axs.plot(data_L_N['L'][ind1], (data_L_N['min_duration_inverse'][ind1] + data_L_N['min_duration_forward'][ind1]) / 2.0, color='black', ls='solid', marker=styles[0])
    p2 = axs.plot(data_L_N['L'][ind2], (data_L_N['min_duration_forward'][ind2] + data_L_N['min_duration_inverse'][ind2]), color='blue', ls='dashed', marker=styles[0])
    lg = axs.legend([p1[0], p2[0], p0[0]], ['Multi resolution', 'Full resolution', 'L$^3$ scaling'], loc='upper left', fontsize=fontsize)
    lg.draw_frame(False)
    axs.set_yscale('log')
    axs.set_xscale('log')
    axs.set_ylim([2e-5, 1e5])
    axs.set_xlim([Lminplot, data_L_N['L'].max() * 1.5])
    axs.set_xlabel('L', fontsize=fontsize)
    axs.set_ylabel('Average duration [s]', fontsize=fontsize)
    # axs.set_title(outname, fontsize=fontsize)
    axs.set_xticks(Lplot)
    axs.set_xticklabels(Lplot)

    fig.tight_layout()
    fig.savefig('s2let_timing_' + suffix + '.pdf', dpi=200)

    fig, axs = plt.subplots(1, 1, figsize=(6, 3.0))

    Lbis = data_L_N['L'][11:]
    p7 = axs.plot(L, 2e-16 * L, color='red')
    p5 = axs.plot(data_L_N['L'][ind1], data_L_N['avg_error'][ind1], color='black', ls='solid', marker=styles[0])
    p6 = axs.plot(data_L_N['L'][ind2], data_L_N['avg_error'][ind2], color='blue', ls='dashed', marker=styles[0])
    lg = axs.legend([p5[0], p6[0], p7[0]], ['Multi resolution', 'Full resolution', 'L scaling'], loc='upper left', fontsize=fontsize)
    lg.draw_frame(False)
    axs.set_yscale('log')
    axs.set_xscale('log')
    axs.set_xlim([Lminplot, data_L_N['L'].max() * 1.5])
    axs.set_xlabel('L', fontsize=fontsize)
    axs.set_ylabel('Maximum error', fontsize=fontsize)
    # axs.set_title(outname, fontsize=fontsize)
    axs.set_ylim([1e-15, 1e-11])
    axs.set_xticks(Lplot)
    axs.set_xticklabels(Lplot)
    
    fig.tight_layout()
    fig.savefig('s2let_error_' + suffix + '.pdf', dpi=200)
    # plt.show()
