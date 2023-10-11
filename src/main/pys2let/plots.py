import numpy as np
import os


def read_matfile(wavelet_type):
    from scipy import io as sio

    matfile = os.path.join(
        os.environ["S2LET"], "src", "main", "matlab", "kappas_" + wavelet_type
    )
    mat_contents = sio.loadmat(matfile)
    kappa = np.ascontiguousarray(mat_contents["kappa_" + wavelet_type])
    kappa0 = np.ascontiguousarray(mat_contents["kappa0_" + wavelet_type].T)
    return kappa, kappa0


def plot_wavs(B, L, J_min, J_max):
    from scipy.interpolate import pchip
    from matplotlib import pyplot as plt
    from pys2let import pys2let_j_max
    import pyssht as ssht

    kappa_spline, kappa0_spline = read_matfile("spline")
    kappa_s2dw, kappa0_s2dw = read_matfile("s2dw")
    kappa_need, kappa0_need = read_matfile("need")

    J = pys2let_j_max(B, L, J_min)

    nx = 1
    ny = 3

    step = 0.01
    xi = np.arange(0, L - 1 + step, step)
    x = np.arange(L)

    plt.figure()
    yi = pchip(x, kappa0_spline)
    plt.semilogx(xi, yi(xi), "-.r")
    yi = pchip(x, kappa0_s2dw)
    plt.semilogx(xi, yi(xi), "-k")
    yi = pchip(x, kappa0_need)
    plt.semilogx(xi, yi(xi), "--b")
    # original loop don't include legend label
    for j in range(J_min, J):
        yi = pchip(x, kappa_spline[j])
        plt.semilogx(xi, yi(xi), "-.r")
        yi = pchip(x, kappa_s2dw[j])
        plt.semilogx(xi, yi(xi), "-k")
        yi = pchip(x, kappa_need[j])
        plt.semilogx(xi, yi(xi), "--b")
    # avoid repeated legend lebels
    yi = pchip(x, kappa_spline[J])
    plt.semilogx(xi, yi(xi), "-.r", label="B-Spline")
    yi = pchip(x, kappa_s2dw[J])
    plt.semilogx(xi, yi(xi), "-k", label="SD")
    yi = pchip(x, kappa_need[J])
    plt.semilogx(xi, yi(xi), "--b", label="Needlet")
    plt.axis([1, L, -0.05, 1.15])
    ticks = 2 ** np.arange(0, J + 3)
    plt.xticks(ticks, ticks)
    plt.xlabel(r"$\ell$")
    plt.legend()

    thetas, phis = ssht.sample_positions(L)
    plt.figure()

    plt.subplot(nx, ny, 1)
    flm = np.zeros(L * L, dtype=complex)
    for l in range(L):
        flm[l * l + l] = kappa0_spline[l]
    f = ssht.inverse(flm, L, Reality=True)
    plt.plot(thetas, f[:, 0], "-.r")
    mx = 1.1 * np.max(f[:, 0])
    plt.axis([0, 2, -mx / 8, mx])
    flm = np.zeros(L * L, dtype=complex)
    for l in range(L):
        flm[l * l + l] = kappa0_s2dw[l]
    f = ssht.inverse(flm, L, Reality=True)
    plt.plot(thetas, f[:, 0], "-k")

    flm = np.zeros(L * L, dtype=complex)
    for l in range(L):
        flm[l * l + l] = kappa0_need[l]
    f = ssht.inverse(flm, L, Reality=True)
    plt.plot(thetas, f[:, 0], "--b")
    plt.xlabel(r"$\theta$")

    for j in range(J_min, J_max + 1):
        plt.subplot(nx, ny, j - J_min + 2)
        flm = np.zeros(L * L, dtype=complex)
        for l in range(L):
            flm[l * l + l] = kappa_spline[j, l]
        f = ssht.inverse(flm, L, Reality=True)
        plt.plot(thetas, f[:, 0], "-.r")
        mx = 1.1 * np.max(f[:, 0])
        plt.axis([0, 2, -mx / 7, mx])
        flm = np.zeros(L * L, dtype=complex)
        for l in range(L):
            flm[l * l + l] = kappa_s2dw[j, l]
        f = ssht.inverse(flm, L, Reality=True)
        plt.plot(thetas, f[:, 0], "-k")
        flm = np.zeros(L * L, dtype=complex)
        for l in range(L):
            flm[l * l + l] = kappa_need[j, l]
        f = ssht.inverse(flm, L, Reality=True)
        plt.plot(thetas, f[:, 0], "--b")
        plt.xlabel(r"$\theta$")

    plt.tight_layout()
    plt.show()


if __name__ == "__main__":
    B = 3
    L = 128
    J_min = 2
    J_max = 3
    plot_wavs(B, L, J_min, J_max)
