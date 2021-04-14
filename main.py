""" main """
from pricing import BlackScholesMerton, GeometricBrownianMotion, Heston
from configuration import ConfigurationBuilder

configuration_obj = ConfigurationBuilder(

    # Common Parameters
    spot=1000.0,
    strike=1000.0,
    sigma=0.20,
    maturity=365,
    risk_free_rate=0,
    dividend_yield=0,

    # Monte carlo Parameters
    simulation=20000,
    steps=365,

    # Heston Parameters
    lt_sigma=0.19,
    rate_reversion=0.1,
    sigma_sigma=1.10,
    correlation=-0.7
)
configuration_obj.strike = 1000
option_1 = BlackScholesMerton(configuration_obj)

configuration_obj.strike = 1100
option_2 = BlackScholesMerton(configuration_obj)

strategy = option_1 * 2

print(option_1.call_price())
print(option_2.call_price())
print(strategy.call_price())
print(strategy.call_delta())






"""
gmb_pricing = GeometricBrownianMotion(configuration_obj)
gmb_pricing.run_simulation()
print(gmb_pricing.call_up_out(barrier=4818))
"""
