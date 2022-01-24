from functools import partial

import numpy as np
from pytest import approx, fixture, mark

from pys2let import (
    analysis_lm2lmn,
    analysis_adjoint_lmn2lm,
    synthesis_lmn2lm,
    synthesis_adjoint_lm2lmn,
    pys2let_j_max,
)


@fixture
def rng(request):
    return np.random.default_rng(getattr(request.config.option, "randomly_seed", None))


def random_lms(rng, L):
    return rng.uniform(size=(L * L, 2)) @ [1, 1j]


def random_wavlet_lms(rng, L, nwvlts):
    lms = [random_lms(rng, L) for _ in range(nwvlts + 1)]
    return lms[0], np.concatenate(lms[1:])


@mark.parametrize(
    "px2wav,wav2px",
    [
        (
            partial(analysis_lm2lmn, spin=0, upsample=1, N=1),
            partial(analysis_adjoint_lmn2lm, spin=0, upsample=1, N=1),
        ),
        (
            partial(analysis_lm2lmn, spin=2, upsample=1, N=1),
            partial(analysis_adjoint_lmn2lm, spin=2, upsample=1, N=1),
        ),
        (
            partial(synthesis_adjoint_lm2lmn, spin=0, upsample=1, N=1),
            partial(synthesis_lmn2lm, spin=0, upsample=1, N=1),
        ),
        (
            partial(synthesis_adjoint_lm2lmn, spin=2, upsample=1, N=1),
            partial(synthesis_lmn2lm, spin=2, upsample=1, N=1),
        ),
    ],
)
def test_axisym_adjoint(
    px2wav, wav2px, rng: np.random.Generator, L=10, B=2, J_min=2
):
    nwvlts = pys2let_j_max(B, L, J_min) - J_min + 1

    x = random_lms(rng, L)
    y_scal, y_wav = random_wavlet_lms(rng, L, nwvlts)

    y = wav2px(y_wav, y_scal, B, L, J_min)
    x_wav, x_scal = px2wav(x, B, L, J_min)

    assert y_wav.conj().T @ x_wav + y_scal.conj() @ x_scal == approx(y.conj().T @ x)
