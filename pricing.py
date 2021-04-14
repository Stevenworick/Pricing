""" pricing module """
import numpy as np
import scipy.stats as stats
from configuration import ConfigurationBuilder


class Option(object):
    """ Option object """
    def __init__(
            self,
            call_price: float,
            put_price: float,
            call_delta: float,
            put_delta: float,
            gamma: float,
            vega: float,
            call_rho: float,
            put_rho: float,
            call_theta: float,
            put_theta: float
    ):
        self._call_price = call_price
        self._put_price = put_price
        self._call_delta = call_delta
        self._put_delta = put_delta
        self._gamma = gamma
        self._vega = vega
        self._call_rho = call_rho
        self._put_rho = put_rho
        self._call_theta = call_theta
        self._put_theta = put_theta

    def call_price(self):
        """ Docstring """
        return self._call_price

    def put_price(self):
        """ Docstring """
        return self._put_price

    def call_delta(self):
        """ Docstring """
        return self._call_delta

    def put_delta(self):
        """ Docstring """
        return self._put_delta

    def gamma(self):
        """ Docstring """
        return self._gamma

    def vega(self):
        """ Docstring """
        return self._vega

    def call_rho(self):
        """ Docstring """
        return self._call_rho

    def put_rho(self):
        """ Docstring """
        return self._put_rho

    def call_theta(self):
        """ Docstring """
        return self._call_theta

    def put_theta(self):
        """ Docstring """
        return self._put_theta


class BlackScholesMerton(object):
    """ Black Scholes Merton Pricing """
    PERIODS_PER_YEAR = 365

    def __init__(self, configuration: ConfigurationBuilder):
        self.configuration = configuration
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

    def __add__(self, other):

        return Option(
            self.call_price() + other.call_price(),
            self.put_price() + other.put_price(),
            self.call_delta() + other.call_delta(),
            self.put_delta() + other.put_delta(),
            self.gamma() + other.gamma(),
            self.vega() + other.vega(),
            self.call_rho() + other.call_rho(),
            self.put_rho() + other.put_rho(),
            self.call_theta() + other.call_theta(),
            self.put_theta() + other.put_theta()
        )


class GeometricBrownianMotion(BlackScholesMerton):
    """ Docstring """
    def __init__(self, configuration: ConfigurationBuilder):
        super().__init__(configuration)
        self.simulation = configuration.simulation
        self.steps = configuration.steps
        self.st_paths = np.zeros((self.simulation, self.steps))
        self.prices_at_maturity = None

    def run_simulation(self):
        """ Docstring """
        for i in range(self.simulation):
            self.st_paths[i][0] = self.s
            for j in range(1, self.steps):
                self.st_paths[i][j] = self.st_paths[i][j-1] * np.exp(
                    (self.r - 0.5 * self.v**2) * 1 / self.PERIODS_PER_YEAR + self.v * np.sqrt(
                        1 / self.PERIODS_PER_YEAR) * np.random.normal(0, 1))
        self.prices_at_maturity = [self.st_paths[i][-1] for i in range(self.simulation)]

    def call_price(self):
        """ Docstring """
        payoffs = [max(S - self.k, 0) for S in self.prices_at_maturity]
        payoff = np.mean(payoffs)
        return payoff * np.exp(-self.r * self.t)

    def put_price(self):
        """ Docstring """
        payoffs = [max(self.k - S, 0) for S in self.prices_at_maturity]
        payoff = np.mean(payoffs)
        return payoff * np.exp(-self.r * self.t)

    def call_up_out(self, barrier):
        """ Docstring """
        payoffs = [(S < barrier) * max(S - self.k, 0) for S in self.prices_at_maturity]
        payoff = np.mean(payoffs)
        return payoff * np.exp(-self.r * self.t)

    def put_down_out(self, barrier):
        """ Docstring """
        payoffs = [(S > barrier) * max(self.k - S, 0) for S in self.prices_at_maturity]
        payoff = np.mean(payoffs)
        return payoff * np.exp(-self.r * self.t)


class Heston(GeometricBrownianMotion):
    """ Stochastic Volatility Pricing """

    def __init__(self, configuration: ConfigurationBuilder):
        super().__init__(configuration)
        self._lt_v = configuration.lt_sigma
        self._rr = configuration.rate_reversion
        self._vv = configuration.sigma_sigma
        self._corr = configuration.correlation
        
    def run_simulation(self):
        """ Docstring """
        delta_t = self.t / float(self.steps)
        for i in range(self.simulation):
            vt = self.v
            st = self.s
            self.st_paths[i][0] = self.s
            for j in range(1, self.steps):
                w1 = np.random.normal(0, 1)
                w2 = self._corr * w1 + np.sqrt(1 - self._corr ** 2) * np.random.normal(0, 1)

                vt = (np.sqrt(max(0, vt)) + 0.5 * self._vv * np.sqrt(delta_t) * w1) ** 2 - self._rr * (
                            vt - self._lt_v) * delta_t - 0.25 * self._vv ** 2 * delta_t

                st = st * np.exp((self.r - self.q - 0.5 * vt) * delta_t + np.sqrt(max(0, vt) * delta_t) * w2)
                self.st_paths[i][j] = st
        self.prices_at_maturity = [self.st_paths[i][-1] for i in range(self.simulation)]
