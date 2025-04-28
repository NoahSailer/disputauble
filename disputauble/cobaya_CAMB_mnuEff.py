import numpy as np

from cobaya import LoggedError
from cobaya.theories.camb import CAMB as CAMBBase
from cobaya.theories.camb.camb import CambTransfers as CambTransfersBase


def difference_dicts(dict1, dict2):
    """
    Compare two dictionaries and return the difference
    """

    assert dict1.keys() == dict2.keys(), "Keys of the dictionaries do not match"
    diff = {}
    for k in dict1.keys():
        if isinstance(dict1[k], (np.ndarray, int, float, complex)) and isinstance(dict2[k], (np.ndarray, int, float, complex)):
            diff[k] = dict1[k] - dict2[k]
        elif isinstance(dict1[k], dict) and isinstance(dict2[k], dict):
            diff[k] = difference_dicts(dict1[k], dict2[k])

    return diff


def multiply_dict_by_factor(dict1, factor):
    """
    Multiply all values in a dictionary by a factor
    """
    out_dict = dict1.copy()
    for k in dict1.keys():
        if isinstance(dict1[k], (np.ndarray, int, float, complex)):
            out_dict[k] *= factor
        elif isinstance(dict1[k], dict):
            out_dict[k] = multiply_dict_by_factor(dict1[k], factor)

    return out_dict


class CobayaCAMB_mnuEff(CAMBBase):

    def calculate(self, state, want_derived=True, **params_values_dict):
        """
        Calculate the theory predictions for the given parameters.

        Parameters
        ----------
        state : dict
            Dictionary to store the results of the calculation.
        want_derived : bool
            Whether to calculate derived parameters.
        params_values_dict : dict
            Dictionary of parameter values.

        Returns
        -------
        None
        """

        transfers = self._camb_transfers.current_state
        if 'results' in transfers:
            self.log.debug("CAMB calculate called with positive neutrino mass")
            super().calculate(state, want_derived=want_derived, **params_values_dict)
        else:
            self.log.debug("CAMB calculate called with negative neutrino mass")
            assert 'results_mnu_abs' in transfers and 'results_mnu_0' in transfers, "CAMB transfers must be calculated before using mnu extrapolation."

            self._camb_transfers.current_state['results'] = transfers['results_mnu_abs']
            mnu_abs_state = {}
            super().calculate(mnu_abs_state, want_derived=want_derived, **params_values_dict)

            self._camb_transfers.current_state['results'] = transfers['results_mnu_0']
            mnu_0_state = {}
            super().calculate(mnu_0_state, want_derived=want_derived, **params_values_dict)

            state.update(difference_dicts(multiply_dict_by_factor(mnu_0_state, 2), mnu_abs_state))

    def must_provide(self, **requirements):
        if 'CAMBdata' in requirements:
            raise LoggedError("CAMBdata object is not supported when extrapolating mnu to negative masses.")

        return super().must_provide(**requirements)

    def get_helper_theories(self):
        """
        Transfer functions are computed separately by camb.transfers, then this
        class uses the transfer functions to calculate power spectra (using A_s, n_s etc).
        """
        self._camb_transfers = CobayaCAMB_transfer_mnuExtrap(self, 'camb.transfers_mnuExtrap',
                                                             dict(stop_at_error=self.stop_at_error),
                                                             timing=self.timer)
        setattr(self._camb_transfers, "requires", self._transfer_requires)
        return {'camb.transfers_mnuExtrap': self._camb_transfers}


class CobayaCAMB_transfer_mnuExtrap(CambTransfersBase):

    def get_can_support_params(self):
        """
        Get the parameters that this theory can support.

        Returns
        -------
        set
            Set of parameters that this theory can support.
        """
        return super().get_can_support_params() | {'mnu_eff'}

    def calculate(self, state, want_derived=True, **params_values_dict):

        if params_values_dict.get('mnu_eff', 0.0) < 0.0:
            self.log.debug("CAMB transfers calculate called with negative neutrino mass")
            mnu_abs_param_values_dict = params_values_dict.copy()
            mnu_0_param_values_dict = params_values_dict.copy()
            if 'mnu_eff' in params_values_dict:
                mnu_abs_param_values_dict['mnu'] = abs(params_values_dict['mnu_eff'])
                del mnu_abs_param_values_dict['mnu_eff']
                self.set_massless_neutrino(mnu_0_param_values_dict)
                del mnu_0_param_values_dict['mnu_eff']

            mnu_0_state = state.copy()
            super().calculate(mnu_0_state, want_derived=want_derived, **mnu_0_param_values_dict)

            mnu_abs_state = state.copy()
            super().calculate(mnu_abs_state, want_derived=want_derived, **mnu_abs_param_values_dict)

            state.update({'results_mnu_0': mnu_0_state['results'], 'results_mnu_abs': mnu_abs_state['results']})
        else:
            if 'mnu_eff' in params_values_dict:
                params_values_dict = params_values_dict.copy()
                params_values_dict['mnu'] = params_values_dict['mnu_eff']
                del params_values_dict['mnu_eff']

            if params_values_dict['mnu'] == 0:
                self.log.debug("CAMB transfers calculate called with massless neutrino")
                self.set_massless_neutrino(params_values_dict)
            else:
                self.log.debug("CAMB transfers calculate called with positive neutrino mass")

            super().calculate(state, want_derived=want_derived, **params_values_dict)

    @staticmethod
    def set_massless_neutrino(param_dict):
        param_dict['mnu'] = 1e-10
        #param_dict['num_nu_massive'] = 0
        #param_dict['nu_mass_eigenstates'] = 0
        #param_dict['nu_mass_degeneracies'] = []
        #param_dict['nu_mass_fractions'] = []
        #param_dict['nu_mass_numbers'] = []
        #param_dict['num_massive_neutrinos'] = 0
