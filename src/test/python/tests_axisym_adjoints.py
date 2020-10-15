from functools import partial

import numpy as np
from pytest import approx, fixture, mark

from pys2let import (
    analysis_adjoint_axisym_wav_mw,
    analysis_adjoint_wav2px,
    analysis_axisym_wav_mw,
    analysis_px2wav,
    synthesis_adjoint_axisym_wav_mw,
    synthesis_adjoint_px2wav,
    synthesis_axisym_wav_mw,
    synthesis_wav2px,
)


@fixture
def rng(request):
    return np.random.default_rng(getattr(request.config.option, "randomly_seed", None))


def random_mw_map(rng, L, spin):
    from pys2let import alm2map_mw

    return alm2map_mw(rng.uniform(size=(L * L, 2)) @ [1, 1j], L, spin)


def random_wavlet_maps(rng, L, spin, nwvlts):
    maps = [random_mw_map(rng, L, spin) for _ in range(nwvlts + 1)]
    return maps[0], np.concatenate(maps[1:])


@mark.parametrize(
    "px2wav,wav2px,spin",
    [
        (analysis_axisym_wav_mw, analysis_adjoint_axisym_wav_mw, 0),
        (synthesis_adjoint_axisym_wav_mw, synthesis_axisym_wav_mw, 0),
        (
            partial(analysis_px2wav, spin=0, upsample=1, N=1),
            partial(analysis_adjoint_wav2px, spin=0, upsample=1, N=1),
            0,
        ),
        (
            partial(analysis_px2wav, spin=2, upsample=1, N=1),
            partial(analysis_adjoint_wav2px, spin=2, upsample=1, N=1),
            2,
        ),
        (
            partial(synthesis_adjoint_px2wav, spin=0, upsample=1, N=1),
            partial(synthesis_wav2px, spin=0, upsample=1, N=1),
            0,
        ),
        (
            partial(synthesis_adjoint_px2wav, spin=2, upsample=1, N=1),
            partial(synthesis_wav2px, spin=2, upsample=1, N=1),
            2,
        ),
    ],
)
def test_axisym_adjoint(
    px2wav, wav2px, spin, rng: np.random.Generator, L=10, B=2, J_min=2
):
    from pys2let import pys2let_j_max

    nwvlts = pys2let_j_max(B, L, J_min) - J_min + 1

    x = random_mw_map(rng, L, spin)
    y_scal, y_wav = random_wavlet_maps(rng, L, spin, nwvlts)

    y = wav2px(y_wav, y_scal, B, L, J_min)
    x_wav, x_scal = px2wav(x, B, L, J_min)

    assert y_wav.conj().T @ x_wav + y_scal.conj() @ x_scal == approx(y.conj().T @ x)
