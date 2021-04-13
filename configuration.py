""" configuration object """


class ConfigurationBuilder(object):
    """ configuration object """
    def __init__(self, spot: float, strike: float, sigma: float, maturity: int, risk_free_rate: float,
                 dividend_yield: float, lt_sigma: float = None, rate_reversion: float = None,
                 sigma_sigma: float = None, correlation: float = None, simulation: int = None, steps: int = None):
        """
        :param spot: current price
        :param strike: exercise price
        :param sigma: current volatility
        :param maturity: number of day until expiration
        :param risk_free_rate: theoretical rate of return of an investment with zero risk
        :param dividend_yield: financial ratio (dividend/price)
        :param lt_sigma: long term average of price volatility
        :param rate_reversion: mean reversion rate for the volatility
        :param sigma_sigma: the volatility of volatility
        :param correlation: correlation between the standard normal random variables W1 and W2
        :param simulation: number of montecarlo simulation
        :param steps: number of steps in each simulation
        """
        self._spot = spot
        self._strike = strike
        self._sigma = sigma
        self._maturity = maturity
        self._risk_free_rate = risk_free_rate
        self._dividend_yield = dividend_yield
        self._lt_sigma = lt_sigma
        self._rate_reversion = rate_reversion
        self._sigma_sigma = sigma_sigma
        self._correlation = correlation
        self._simulation = simulation
        self._steps = steps

    @property
    def spot(self):
        """ getter """
        return self._spot

    @spot.setter
    def spot(self, _spot):
        """ setter """
        self._spot = _spot

    @property
    def strike(self):
        """ getter """
        return self._strike

    @strike.setter
    def strike(self, _strike):
        """ setter """
        self._strike = _strike

    @property
    def sigma(self):
        """ getter """
        return self._sigma

    @sigma.setter
    def sigma(self, _sigma):
        """ setter """
        self._sigma = _sigma

    @property
    def maturity(self):
        """ getter """
        return self._maturity

    @maturity.setter
    def maturity(self, _maturity):
        """ setter """
        self._maturity = _maturity

    @property
    def risk_free_rate(self):
        """ getter """
        return self._risk_free_rate

    @risk_free_rate.setter
    def risk_free_rate(self, _risk_free_rate):
        """ setter """
        self._risk_free_rate = _risk_free_rate

    @property
    def dividend_yield(self):
        """ getter """
        return self._dividend_yield

    @dividend_yield.setter
    def dividend_yield(self, _dividend_yield):
        """ setter """
        self._dividend_yield = _dividend_yield

    @property
    def lt_sigma(self):
        """ getter """
        return self._lt_sigma

    @lt_sigma.setter
    def lt_sigma(self, _lt_sigma):
        """ setter """
        self._lt_sigma = _lt_sigma

    @property
    def rate_reversion(self):
        """ getter """
        return self._rate_reversion

    @rate_reversion.setter
    def rate_reversion(self, _rate_reversion):
        """ setter """
        self._rate_reversion = _rate_reversion

    @property
    def sigma_sigma(self):
        """ getter """
        return self._sigma_sigma

    @sigma_sigma.setter
    def sigma_sigma(self, _sigma_sigma):
        """ setter """
        self._sigma_sigma = _sigma_sigma

    @property
    def correlation(self):
        """ getter """
        return self._correlation

    @correlation.setter
    def correlation(self, _correlation):
        """ setter """
        self._correlation = _correlation

    @property
    def simulation(self):
        """ getter """
        return self._simulation

    @simulation.setter
    def simulation(self, _simulation):
        """ setter """
        self._simulation = _simulation

    @property
    def steps(self):
        """ getter """
        return self._steps

    @steps.setter
    def steps(self, _steps):
        """ setter """
        self._steps = _steps