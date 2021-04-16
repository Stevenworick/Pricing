""" Exotic Option Replication : Call Up & Out """
from pricing import BlackScholesMerton, GeometricBrownianMotion
from configuration import ConfigurationBuilder

configuration = ConfigurationBuilder(
    kind='put',
    spot=100.0,
    sigma=0.20,
    maturity=365,
    risk_free_rate=0,
    dividend_yield=0,
)

# :::: Buy :::: ATM :::: Call
configuration.strike = configuration.spot
option_1 = BlackScholesMerton(configuration)

# :::: Sell :::: OTM :::: Call :::: Barrier @ 110
barrier = 110
configuration.strike = barrier
option_2 = BlackScholesMerton(configuration)

# :::: Buy :::: OTM :::: Call :::: Shift @ 1
shift = 0.001  # int(barrier * 0.03)
configuration.strike = barrier + shift
option_3 = BlackScholesMerton(configuration)

# :::: Aggregate
q = (barrier-configuration.spot)/shift
strategy = option_1 - option_2 * (q + 1) + option_3 * q

print(strategy.price())
print(strategy.delta())
print(strategy.gamma())
print(strategy.vega())
print(strategy.theta())
print(strategy.rho())

"""
Result with shift @ 1
0.636295430446495
0.00786088027331111
-0.0006840904252101709
-0.02052271275630435
0.0008433991543687069
0.0014979259688452373

Result with shift @ 0.001
0.5825870460103033
0.0070055752839834895
-0.0006294851178267891
-0.01888455353491736
0.0007760775425538213
0.0011797048241533048
"""

# :::: Montecarlo for comparison
configuration.strike = configuration.spot
configuration.simulation = 40000
configuration.steps = 365
gbm_pricing = GeometricBrownianMotion(configuration)
gbm_pricing.run_simulation()
print(gbm_pricing.call_up_out(barrier=barrier))

"""
Result with Montecarlo
0.5858203229501271
"""