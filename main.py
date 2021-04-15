""" main """
from pricing import BlackScholesMerton, GeometricBrownianMotion, Heston
from configuration import ConfigurationBuilder

configuration = ConfigurationBuilder(
    kind='call',
    spot=100.0,
    sigma=0.30,
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
shift = int(barrier * 0.03)
configuration.strike = barrier + shift
option_3 = BlackScholesMerton(configuration)

# :::: Aggregate
q = (barrier-configuration.spot)/shift
strategy = option_1 - option_2 * (q + 1) + option_3 * q


print(strategy.price())
print(strategy.delta())

