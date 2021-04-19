""" configuration object """


class OptionConfigurationBuilder(object):
    """ configuration object """
    def __init__(self, kind: str = None, spot: float = None, strike: float = None, sigma: float = None,
                 maturity: int = None, risk_free_rate: float = None, dividend_yield: float = None,
                 simulation: int = None, steps: int = None):
        """
        :param kind: should be 'call' or 'put'
        :param spot: current price
        :param strike: exercise price
        :param sigma: current volatility
        :param maturity: number of day until expiration
        :param risk_free_rate: theoretical rate of return of an investment with zero risk
        :param dividend_yield: financial ratio (dividend/price)
        :param simulation: number of montecarlo simulation
        :param steps: number of steps in each simulation
        """
        self._kind = kind
        self._spot = spot
        self._strike = strike
        self._sigma = sigma
        self._maturity = maturity
        self._risk_free_rate = risk_free_rate
        self._dividend_yield = dividend_yield
        self._simulation = simulation
        self._steps = steps

    @property
    def kind(self):
        """ getter """
        return self._kind

    @kind.setter
    def kind(self, _kind):
        """ setter """
        self._kind = _kind

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
