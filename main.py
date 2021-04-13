""" main """

from configuration import ConfigurationBuilder
from pricing import BSM, Heston


configuration_obj = ConfigurationBuilder(
    spot=100.0,
    strike=100.0,
    sigma=0.25,
    maturity=252,
    risk_free_rate=0,
    dividend_yield=0,
    lt_sigma=0.15,
    rate_reversion=0.1,
    sigma_sigma=0.15,
    correlation=-0.1,
    simulation=10000,
    steps=252
)

vanilla_option = BSM(configuration_obj)
exotic_option = Heston(configuration_obj)
