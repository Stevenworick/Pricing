""" pricing module """
import numpy as np
import scipy.stats as stats
from configuration import ConfigurationBuilder


class BSM(object):
    """ Black Scholes Merton Pricing """
    PERIODS_PER_YEAR = 365

    def __init__(self, configuration: ConfigurationBuilder):
        self._s = configuration.spot
        self._k = configuration.strike
        self._v = configuration.sigma
        self._t = configuration.maturity / self.PERIODS_PER_YEAR
        self._r = configuration.risk_free_rate
        self._q = configuration.dividend_yield
        self.d1 = (np.log(self._s / self._k) + (self._r + self._v ** 2 / 2.0) * self._t) / self._v * np.sqrt(self._t)
        self.d2 = self.d1 - self._v * np.sqrt(self._t)

    def call_price(self):
        """ Docstring """
        price = np.exp(-self._r * self._t) * (self._s * np.exp((self._r - self._q) * self._t) * stats.norm.cdf(
            self.d1) - self._k * stats.norm.cdf(self.d2))
        return price

    def put_price(self):
        """ Docstring """
        price = np.exp(-self._r * self._t) * (self._k * stats.norm.cdf(-self.d2) - (
                    self._s * np.exp((self._r - self._q) * self._t) * stats.norm.cdf(-self.d1)))

        return price

    def call_delta(self):
        """ Docstring """
        return np.exp(-self._q * self._t) * stats.norm.cdf(self.d1)

    def put_delta(self):
        """ Docstring """
        return np.exp(-self._q * self._t) * stats.norm.cdf(self.d1) - 1

    def gamma(self):
        """ Docstring """
        return stats.norm.pdf(self.d1) * np.exp(-self._q * self._t) / (self._s * self._v * np.sqrt(self._t))

    def vega(self):
        """ Docstring """
        return 0.01 * (self._s * np.sqrt(self._t) * stats.norm.pdf(self.d1) * np.exp(-self._q * self._t))

    def call_rho(self):
        """ Docstring """
        return 0.01 * (self._k * self._t * (np.exp(-self._r * self._t)) * stats.norm.cdf(self.d2))

    def put_rho(self):
        """ Docstring """
        return 0.01 * (-self._k * self._t * (np.exp(-self._r * self._t)) * stats.norm.cdf(-self.d2))

    def call_theta(self):
        """ Docstring """
        theta = -self._s * stats.norm.pdf(self.d1) * self._v * np.exp(-self._q * self._t) / (2 * np.sqrt(self._t)) \
                + self._q * self._s * stats.norm.cdf(self.d1) * np.exp(-self._q * self._t) \
                - self._r * self._k * np.exp(-self._r * self._t) * stats.norm.cdf(self.d2)
        return 1/self.PERIODS_PER_YEAR * theta

    def put_theta(self):
        """ Docstring """
        theta = -self._s * stats.norm.pdf(self.d1) * self._v * np.exp(-self._q * self._t) / (2 * np.sqrt(self._t)) \
                + self._q * self._s * stats.norm.cdf(-self.d1) * np.exp(-self._q * self._t) \
                - self._r * self._k * np.exp(-self._r * self._t) * stats.norm.cdf(-self.d2)
        return 1/self.PERIODS_PER_YEAR * theta


class Heston(BSM):
    """ Stochastic Volatility Pricing """
    PERIODS_PER_YEAR = 365

    def __init__(self, configuration: ConfigurationBuilder):
        super().__init__(configuration)
        self._lt_v = configuration.lt_sigma
        self._rr = configuration.rate_reversion
        self._vv = configuration.sigma_sigma
        self._corr = configuration.correlation
        self._simulation = configuration.simulation
        self._steps = configuration.steps
        self._st_paths = np.zeros((self._simulation, self._steps))
        self.run_simulations()

    @property
    def st_paths(self):
        """ Docstring """
        return self._st_paths

    def run_simulations(self):
        """ Docstring """
        delta_t = self._t / float(self._steps)
        for i in range(self._simulation):
            vt = self._v
            st = self._s
            self._st_paths[i][0] = self._s
            for j in range(1, self._steps):
                w1 = np.random.normal(0, 1)
                w2 = self._corr * w1 + np.sqrt(1 - self._corr ** 2) * np.random.normal(0, 1)
                vt = (np.sqrt(max(0, vt)) + 0.5 * self._vv * np.sqrt(delta_t) * w1) ** 2 - self._rr * (
                            vt - self._lt_v) * delta_t - 0.25 * self._vv ** 2 * delta_t

                st = st * np.exp((self._r - 0.5 * vt) * delta_t + np.sqrt(max(0, vt) * delta_t) * w2)
                self._st_paths[i][j] = st

    def call_price(self):
        """ Docstring """
        return 0
