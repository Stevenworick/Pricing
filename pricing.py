""" pricing module """
import numpy as np
import scipy.stats as stats
from configuration import ConfigurationBuilder


class Option(object):
    """ Option object """
    def __init__(self, price: float = None, delta: float = None, gamma: float = None, vega: float = None,
                 theta: float = None, rho: float = None):
        self._price = price
        self._delta = delta
        self._gamma = gamma
        self._vega = vega
        self._theta = theta
        self._rho = rho

    def price(self):
        """ Docstring """
        return self._price

    def delta(self):
        """ Docstring """
        return self._delta

    def gamma(self):
        """ Docstring """
        return self._gamma

    def vega(self):
        """ Docstring """
        return self._vega

    def theta(self):
        """ Docstring """
        return self._theta

    def rho(self):
        """ Docstring """
        return self._rho

    def __add__(self, other):
        """ Docstring """
        return Option(
            self.price() + other.price(),
            self.delta() + other.delta(),
            self.gamma() + other.gamma(),
            self.vega() + other.vega(),
            self.theta() + other.theta(),
            self.rho() + other.rho()
        )

    def __sub__(self, other):
        """ Docstring """
        return Option(
            self.price() - other.price(),
            self.delta() - other.delta(),
            self.gamma() - other.gamma(),
            self.vega() - other.vega(),
            self.theta() - other.theta(),
            self.rho() - other.rho()
        )

    def __mul__(self, other):
        """ Docstring """
        return Option(
            self.price() * other,
            self.delta() * other,
            self.gamma() * other,
            self.vega() * other,
            self.theta() * other,
            self.rho() * other
        )


class BlackScholesMerton(Option):
    """ Black Scholes Merton Pricing """
    PERIODS_PER_YEAR = 365

    def __init__(self, configuration: ConfigurationBuilder):
        super(BlackScholesMerton, self).__init__()
        self.kind = configuration.kind
        self.s = configuration.spot
        self.k = configuration.strike
        self.v = configuration.sigma
        self.t = configuration.maturity / self.PERIODS_PER_YEAR
        self.r = configuration.risk_free_rate
        self.q = configuration.dividend_yield
        self._d1 = (np.log(self.s / self.k) + (self.r + self.v ** 2 / 2.0) * self.t) / self.v * np.sqrt(self.t)
        self._d2 = self._d1 - self.v * np.sqrt(self.t)

    def price(self):
        """ Docstring """
        if self.kind == 'call':
            price = np.exp(-self.r * self.t) * (self.s * np.exp((self.r - self.q) * self.t) * stats.norm.cdf(
                self._d1) - self.k * stats.norm.cdf(self._d2))
        elif self.kind == 'put':
            price = np.exp(-self.r * self.t) * (self.k * stats.norm.cdf(-self._d2) - (
                    self.s * np.exp((self.r - self.q) * self.t) * stats.norm.cdf(-self._d1)))
        else:
            raise ValueError
        return price

    def delta(self):
        """ Docstring """
        if self.kind == 'call':
            delta = np.exp(-self.q * self.t) * stats.norm.cdf(self._d1)
        elif self.kind == "put":
            delta = np.exp(-self.q * self.t) * stats.norm.cdf(self._d1) - 1
        else:
            raise ValueError
        return delta

    def gamma(self):
        """ Docstring """
        return stats.norm.pdf(self._d1) * np.exp(-self.q * self.t) / (self.s * self.v * np.sqrt(self.t))

    def vega(self):
        """ Docstring """
        return 0.01 * (self.s * np.sqrt(self.t) * stats.norm.pdf(self._d1) * np.exp(-self.q * self.t))

    def rho(self):
        """ Docstring """
        if self.kind == 'call':
            rho = 0.01 * (self.k * self.t * (np.exp(-self.r * self.t)) * stats.norm.cdf(self._d2))
        elif self.kind == "put":
            rho = 0.01 * (-self.k * self.t * (np.exp(-self.r * self.t)) * stats.norm.cdf(-self._d2))
        else:
            raise ValueError
        return rho

    def theta(self):
        """ Docstring """
        if self.kind == 'call':
            theta = -self.s * stats.norm.pdf(self._d1) * self.v * np.exp(-self.q * self.t) / (2 * np.sqrt(self.t)) \
                    + self.q * self.s * stats.norm.cdf(self._d1) * np.exp(-self.q * self.t) \
                    - self.r * self.k * np.exp(-self.r * self.t) * stats.norm.cdf(self._d2)
        elif self.kind == "put":
            theta = -self.s * stats.norm.pdf(self._d1) * self.v * np.exp(-self.q * self.t) / (2 * np.sqrt(self.t)) \
                    + self.q * self.s * stats.norm.cdf(-self._d1) * np.exp(-self.q * self.t) \
                    - self.r * self.k * np.exp(-self.r * self.t) * stats.norm.cdf(-self._d2)
        else:
            raise ValueError
        return 1/self.PERIODS_PER_YEAR * theta


class GeometricBrownianMotion(BlackScholesMerton):
    """ Docstring """
    def __init__(self, configuration: ConfigurationBuilder):
        super(GeometricBrownianMotion, self).__init__(configuration)
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

    def price(self):
        """ Docstring """
        if self.kind == 'call':
            payoffs = [max(S - self.k, 0) for S in self.prices_at_maturity]
        elif self.kind == 'put':
            payoffs = [max(self.k - S, 0) for S in self.prices_at_maturity]
        else:
            raise ValueError
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
