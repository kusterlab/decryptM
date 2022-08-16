# filters.py
#
# Resource file for regulated curve filters
#
# Florian P. Bayer - Aug. 2022
#
# version 1.0
#

# filter for down curves from ms2 data
ms2_down = {
    # Minimal logistic R2 fit threshold
    'R2 MINIMUM': 0.7,
    # The Maximal Root Mean Squared Error of the Curve fit to the observations
    'RMSE MAXIMUM': 0.25,
    # Include curves if curve top fit is inside range (min, max)
    'TOP RANGE': (0.7, 1/0.7),
    # Exclude curves if bottom inside range (min, max)
    'BOTTOM MAXIMUM': 0.6,
    # Maximum effect
    'EFFECT SIZE MAXIMUM': -0.4,
    # Exclude curves if slope is steeper than slope max
    'SLOPE MAXIMUM': 6,
    # Exclude curves if high dose effect higher than max high dose
    'HIGH DOSE MAXIMUM': 0.6,
    # Exclude curves if lowest dose effect outside range (min, max)
    'FIRST DOSE RANGE': (0.6, 1.66),
}

# filter for down curves from ms3 data
ms3_down = {
    # Minimal logistic R2 fit threshold
    'R2 MINIMUM': 0.7,
    # The Maximal Root Mean Squared Error of the Curve fit to the observations
    'RMSE MAXIMUM': 0.25,
    # Include curves if curve top fit is inside range (min, max)
    'TOP RANGE': (0.7, 1/0.7),
    # Exclude curves if bottom inside range (min, max)
    'BOTTOM MAXIMUM': 0.5,
    # Maximum effect required
    'EFFECT SIZE MAXIMUM': -0.5,
    # Exclude curves if slope is steeper than slope max
    'SLOPE MAXIMUM': 6,
    # Exclude curves if high dose effect higher than max high dose
    'HIGH DOSE MAXIMUM': 0.5,
    # Exclude curves if lowest dose effect outside range (min, max)
    'FIRST DOSE RANGE': (0.5, 2),
}

# filter for up curves from ms2 data
# First dose and curve top as a function of the effect size
ms2_up = {
    # Minimal logistic R2 fit threshold
    'R2 MINIMUM': 0.7,
    # The Maximal Relative Root Mean Squared Error of the Curve fit to the observations
    'RMSE MAXIMUM': lambda mean_ratio: mean_ratio * 0.25,
    # Exclude curves if  bigger than dynamic cut of top value
    'TOP MAXIMUM': lambda effect: effect * 0.05 + 1.35,
    # Exclude curves if bottom is smaller than minimum value
    'BOTTOM MINIMUM': 1.66,
    # Exclude curves if effect size is smaller than minimum value
    'EFFECT SIZE MINIMUM': 0.66,
    # Exclude curves if slope is steeper than slope max
    'SLOPE MAXIMUM': 6,
    # Exclude curves if high dose effect inside range (min, max)
    'HIGH DOSE MINIMUM': 1.66,
    # Exclude curves if lowest dose effect higher than dynamic cut of highest dose
    'FIRST DOSE MAXIMUM': lambda effect: effect * 0.05 + 1.35,
}

# filter for up curves from ms3 data
# First dose and curve top as a function of the effect size
ms3_up = {
    # Minimal logistic R2 fit threshold
    'R2 MINIMUM': 0.7,
    # The Maximal Relative Root Mean Squared Error of the Curve fit to the observations
    'RMSE MAXIMUM': lambda mean_ratio: mean_ratio * 0.25,
    # Exclude curves if  bigger than dynamic cut of top value
    'TOP MAXIMUM': lambda effect: effect * 0.05 + 1.35,
    # Exclude curves if bottom is smaller than minimum value
    'BOTTOM MINIMUM': 2,
    # Exclude curves if effect size is smaller than minimum value
    'EFFECT SIZE MINIMUM': 1,
    # Exclude curves if slope is steeper than slope max
    'SLOPE MAXIMUM': 6,
    # Exclude curves if high dose effect inside range (min, max)
    'HIGH DOSE MINIMUM': 2,
    # Exclude curves if lowest dose effect higher than dynamic cut of highest dose
    'FIRST DOSE MAXIMUM': lambda effect: effect * 0.05 + 1.35,
}
