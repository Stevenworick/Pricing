""" main """

from configuration import ConfigurationBuilder
from pricing import BSM, Heston

configuration_obj = ConfigurationBuilder(

    # Common Parameters
    spot=100.0,
    strike=100.0,
    sigma=0.25,
    maturity=252,
    risk_free_rate=0,
    dividend_yield=0,

    # Heston Parameters
    lt_sigma=0.19,
    rate_reversion=0,
    sigma_sigma=1.10,
    correlation=-0.7,
    simulation=1000,
    steps=252
)

bsm_pricing = BSM(configuration_obj)
heston_pricing = Heston(configuration_obj)

print(bsm_pricing.call_price())
print(bsm_pricing.put_price())
print(heston_pricing.call_price())
print(heston_pricing.put_price())
