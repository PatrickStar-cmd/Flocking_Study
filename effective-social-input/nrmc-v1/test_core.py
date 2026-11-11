"""Small, data-free NRMC contract check. Run: python test_core.py."""
import numpy as np

from model import Field
from nrmc import Bank, compensate, candidates


def main():
    model = dict(N=4, M=40, L=1000., dt=.3, beta=1000., speed_scale=.05,
                 global_inhibition=0., connection_exponent=.5,
                 kernel_sigma=.4, edge_weight=.012)
    filt = dict(R=1e-4, q_H=1e-4, q_C=1e-4, q_D=1e-5, tau=10.,
                epsilon=.02, samples=64, gamma=np.arange(0., 1.001, .05))
    field = Field(model)
    state = field.initial(321)
    full = np.ones((field.n, field.n-1), dtype=bool)
    y = field.bearings(state)
    bank = Bank(full.shape, filt, model['dt'])
    bank.update(y, full)
    exact = field.observed_input(state['alpha'], y, full)
    got, _ = compensate('nrmc', bank, y, full, state['u'], state['alpha'],
                        field, filt, [321, 0])
    assert np.array_equal(got, exact)

    fresh = full.copy()
    fresh[:, 1:] = False
    hidden_y = y.copy()
    hidden_y[~fresh] = np.nan
    bank.update(hidden_y, fresh)
    args = ('nrmc', bank, hidden_y, fresh, state['u'], state['alpha'],
            field, filt, [321, 1])
    output, diag = compensate(*args)
    again, again_diag = compensate(*args)
    assert np.array_equal(output, again) and diag == again_diag
    assert np.isfinite(output).all()
    assert diag['selected_risk'] <= diag['mix_risk'] + 1e-12
    assert candidates(np.zeros((field.n, 4, field.m)), filt['gamma'])[0].shape[1] == 81

    leaked = hidden_y.copy()
    leaked[~fresh] = 0.
    try:
        bank.copy().update(leaked, fresh)
    except ValueError:
        pass
    else:
        raise AssertionError('Hidden bearings accepted')
    print('PASS: fresh input, hidden observation contract, deterministic NRMC, 81 candidates')


if __name__ == '__main__':
    main()
