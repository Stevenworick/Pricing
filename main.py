""" main """
from configuration import OptionConfigurationBuilder
from pricing import BlackScholesMerton

if __name__ == '__main__':

    configuration = OptionConfigurationBuilder(
        kind='call',
        spot=100,
        strike=100,
        sigma=0.5,
        maturity=252,
        risk_free_rate=0,
        dividend_yield=0,
        simulation=10000,
        steps=252
    )
    test = BlackScholesMerton(configuration)
    test_ = BlackScholesMerton(configuration)
    print('hash is', hash(test) == hash(test_))