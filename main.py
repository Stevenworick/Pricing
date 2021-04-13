""" main """

from configuration import ConfigurationBuilder
from pricing import BlackScholesMerton, GeometricBrownianMotion, Heston

configuration_obj = ConfigurationBuilder(

    # Common Parameters
    spot=100.0,
    strike=100.0,
    sigma=0.25,
    maturity=252,
    risk_free_rate=0,
    dividend_yield=0,

    # Monte carlo Parameters
    simulation=20000,
    steps=252,

    # Heston Parameters
    lt_sigma=0.19,
    rate_reversion=0.1,
    sigma_sigma=1.10,
    correlation=-0.7
)

bsm_pricing = BlackScholesMerton(configuration_obj)
gbm_pricing = GeometricBrownianMotion(configuration_obj)
heston_pricing = Heston(configuration_obj)

print("call")
print(bsm_pricing.call_price())
print(gbm_pricing.call_price())
print(heston_pricing.call_price())

print("put")
print(bsm_pricing.put_price())
print(gbm_pricing.put_price())
print(heston_pricing.put_price())

"""
In[2]:
9.94764496602258
10.089992268159822
9.96097452176407
put
9.94764496602258
9.925691552994227
9.876303619009215
"""