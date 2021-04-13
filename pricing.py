""" pricing module """
import numpy as np
import scipy.stats as stats
from configuration import ConfigurationBuilder


class BSM(object):
    """ Black Scholes Merton Pricing """
    PERIODS_PER_YEAR = 252

    def __init__(self, configuration: ConfigurationBuilder):
        self.s = configuration.spot
        self.k = configuration.strike
        self.v = configuration.sigma
        self.t = configuration.maturity / self.PERIODS_PER_YEAR
        self.r = configuration.risk_free_rate
        self.q = configuration.dividend_yield
        self._d1 = (np.log(self.s / self.k) + (self.r + self.v ** 2 / 2.0) * self.t) / self.v * np.sqrt(self.t)
        self._d2 = self._d1 - self.v * np.sqrt(self.t)

    def call_price(self):
        """ Docstring """
        price = np.exp(-self.r * self.t) * (self.s * np.exp((self.r - self.q) * self.t) * stats.norm.cdf(
            self._d1) - self.k * stats.norm.cdf(self._d2))
        return price

    def put_price(self):
        """ Docstring """
        price = np.exp(-self.r * self.t) * (self.k * stats.norm.cdf(-self._d2) - (
                    self.s * np.exp((self.r - self.q) * self.t) * stats.norm.cdf(-self._d1)))

        return price

    def call_delta(self):
        """ Docstring """
        return np.exp(-self.q * self.t) * stats.norm.cdf(self._d1)

    def put_delta(self):
        """ Docstring """
        return np.exp(-self.q * self.t) * stats.norm.cdf(self._d1) - 1

    def gamma(self):
        """ Docstring """
        return stats.norm.pdf(self._d1) * np.exp(-self.q * self.t) / (self.s * self.v * np.sqrt(self.t))

    def vega(self):
        """ Docstring """
        return 0.01 * (self.s * np.sqrt(self.t) * stats.norm.pdf(self._d1) * np.exp(-self.q * self.t))

    def call_rho(self):
        """ Docstring """
        return 0.01 * (self.k * self.t * (np.exp(-self.r * self.t)) * stats.norm.cdf(self._d2))

    def put_rho(self):
        """ Docstring """
        return 0.01 * (-self.k * self.t * (np.exp(-self.r * self.t)) * stats.norm.cdf(-self._d2))

    def call_theta(self):
        """ Docstring """
        theta = -self.s * stats.norm.pdf(self._d1) * self.v * np.exp(-self.q * self.t) / (2 * np.sqrt(self.t)) \
                + self.q * self.s * stats.norm.cdf(self._d1) * np.exp(-self.q * self.t) \
                - self.r * self.k * np.exp(-self.r * self.t) * stats.norm.cdf(self._d2)
        return 1/self.PERIODS_PER_YEAR * theta

    def put_theta(self):
        """ Docstring """
        theta = -self.s * stats.norm.pdf(self._d1) * self.v * np.exp(-self.q * self.t) / (2 * np.sqrt(self.t)) \
                + self.q * self.s * stats.norm.cdf(-self._d1) * np.exp(-self.q * self.t) \
                - self.r * self.k * np.exp(-self.r * self.t) * stats.norm.cdf(-self._d2)
        return 1/self.PERIODS_PER_YEAR * theta


class Heston(BSM):
    """ Stochastic Volatility Pricing """

    def __init__(self, configuration: ConfigurationBuilder):
        super().__init__(configuration)
        self._lt_v = configuration.lt_sigma
        self._rr = configuration.rate_reversion
        self._vv = configuration.sigma_sigma
        self._corr = configuration.correlation
        self._simulation = configuration.simulation
        self._steps = configuration.steps
        self._st_paths = np.zeros((self._simulation, self._steps))
        self._prices_at_maturity = None
        self.run_simulations()

    @property
    def st_paths(self):
        """ Docstring """
        return self._st_paths

    @property
    def prices_at_maturity(self):
        """ Docstring """
        return self._prices_at_maturity

    def run_simulations(self):
        """ Docstring """
        delta_t = self.t / float(self._steps)
        for i in range(self._simulation):
            vt = self.v
            st = self.s
            self._st_paths[i][0] = self.s
            for j in range(1, self._steps):
                w1 = np.random.normal(0, 1)
                w2 = self._corr * w1 + np.sqrt(1 - self._corr ** 2) * np.random.normal(0, 1)

                vt = (np.sqrt(max(0, vt)) + 0.5 * self._vv * np.sqrt(delta_t) * w1) ** 2 - self._rr * (
                            vt - self._lt_v) * delta_t - 0.25 * self._vv ** 2 * delta_t

                st = st * np.exp((self.r - self.q - 0.5 * vt) * delta_t + np.sqrt(max(0, vt) * delta_t) * w2)
                self._st_paths[i][j] = st
        self._prices_at_maturity = [self._st_paths[i][-1] for i in range(self._simulation)]

    def call_price(self):
        """ Docstring """
        payoffs = [max(S - self.k, 0) for S in self._prices_at_maturity]
        payoff = np.mean(payoffs)
        return payoff * np.exp(-self.r * self.t)

    def put_price(self):
        """ Docstring """
        payoffs = [max(self.k - S, 0) for S in self._prices_at_maturity]
        payoff = np.mean(payoffs)
        return payoff * np.exp(-self.r * self.t)
