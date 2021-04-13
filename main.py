""" main """
from pricing import BlackScholesMerton, GeometricBrownianMotion, Heston
from configuration import ConfigurationBuilder

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
gbm_pricing.run_simulation()

heston_pricing = Heston(configuration_obj)
heston_pricing.run_simulation()

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
call
9.94764496602258
10.090981528256577
2.3511986201121142
put
9.94764496602258
9.966650183954462
29.61329599783158
"""