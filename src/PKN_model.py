import math
import pandas as pd
from DeclineFunction import small_g, exponentAlfa

PI = 3.1416

def TSO_extended_time(mult, t_tso, eff_tso, n_prime_rheo, C_L, spurt_gp100ft):
    """
    Used to calculate extended time for tip screen out.

    Params
    -------
    mult: float
        TSO multiplier, obtained from inverse function of TSO_width_multiplier
    t_tso: float
        Pumping time (also known as tp)
    eff_tso: float
        Pumping efficiency (vol. of fracture/total vol.)
    C_L: float
        Fluid loss coefficient
    spurt_gp100ft: float
        Spurt loss

    Returns
    -------
    tso_extend_time: float
        Time extended during tip screen out.
    """
    add_tp_test = t_tso * 0.5
    spurt = spurt_gp100ft

    for i in range(50): # use finite loop to avoid unconvergence
        test_mult = TSO_width_multiplier(
            add_tp_test, 
            t_tso, 
            eff_tso, 
            n_prime_rheo, 
            C_L, 
            spurt
        )
        dif = (test_mult - mult) / mult

        if abs(dif) <= 0.01:
            break
        else:
            if dif < 0:
                add_tp_test = add_tp_test * (1 + abs(dif))
            else:
                add_tp_test = add_tp_test * (1 - abs(dif))

    tso_extend_time = add_tp_test
    
    return tso_extend_time

def TSO_width_multiplier(add_tp, t_tso, eff_tso, n_prime_rheo, C_L, spurt_gp100ft):
    """
    Calculate needed multiplier for tip screen out width.

    Params
    -------
    add_tp: float
        Additional pump time
    t_tso: float
        Pumping time (also known as tp)
    eff_tso: float
        Pumping efficiency (vol. of fracture/total vol.)
    C_L: float
        Fluid loss coefficient
    spurt_gp100ft: float
        Spurt loss

    Returns
    -------
    TSO_width_mult: float
        Multiplier for tip screen out width
    """
    spurt = spurt_gp100ft * 0.001336806    # gal/100ft2 -> ft

    g0 = small_g(0, n_prime_rheo, eff_tso)
    dtD = add_tp / t_tso
    g_dtD = small_g(dtD, n_prime_rheo, eff_tso)

    kappa = 1 + spurt / (g0 * C_L * math.sqrt(t_tso))

    TSO_width_mult = 1 + (dtD / eff_tso) - (1 - eff_tso) / (kappa * eff_tso) * (g_dtD / g0 - 1)

    return TSO_width_mult

def modelPKN_tp(frac_length, q, frac_height, fl_height, w_avg_in, cl, spurt_gp100ft):
    """
    Calculates the pumping time (tp)

    Params
    -------
    frac_length: float
        Fracture design length
    q: float
        Pump speed
    frac_height: float
        Fracture height
    fl_height: float
        Fluid loss height
    w_avg_in: float
        Average fracture width in inch
    cl: float
        Fluid loss coefficient
    spurt_gp100ft:
        Spurt loss

    Returns
    -------
    PKN_tp: float
        Pumping time
    """
    spurt = spurt_gp100ft * 0.001336806    # gal/100ft2 -> ft
    w_avg_ft = w_avg_in / 12               # in -> ft
    r_p = fl_height / frac_height
    A_f = 2 * frac_height * frac_length

    go = 1.5
    vol_lost = 2 * r_p * A_f * go * cl

    # solution of equation of second order (y= ax2+bx+c) where x=sqr(t)
    a = q
    b = -3 * A_f * r_p * cl
    c = -A_f * (w_avg_ft + 2 * r_p * spurt)

    sqr_t = (-b + math.sqrt(b**2 - 4 * a * c)) / (2 * a)
    
    return sqr_t**2

def frac_length(Vi, Q, frac_ht, fl_ht, plain_strain, n, K, cl, spurt_gp100ft):
    """
    Find fracture length based on PKN model

    Params
    -------
    Vi: float
        Volume of injected fluid
    frac_ht: float
        Fracture height
    fl_ht: float
        Fluid loss height
    plain_strain: float
        Plaint strain modulus (in lb-s^n/ft2)
    n: float
        Power Law rheology constants
    K: float
        Another Power Law rheology constants
    cl: float
        Fluid loss coefficient
    spurt_gp100ft: float
        Spurt loss

    Returns
    -------
    PKN_length: float
        Length of fracture calculated with PKN model
    """
    w_avg_ft = w_avg_in / 12
    spurt_ft = spurt_gp100ft * 0.001336806
    q_ft3 = Q * 5.614583

    length_test = 500 # in ft
    eff_tmp = 0.5
    alfa = exponentAlfa(n, eff_tmp)
    for _ in range(1,50):
        w_max = 3 * ((K * length_test) * (frac_ht / plain_strain) * (q_ft3 / (frac_ht * 60))**n) ** (1 / (2 * n+2))
        beta_p = (n + 2) / (n + 3)
        w_avg_in = 12 * PI / 4 * beta_p * w_max
        A_f = 2 * frac_ht * length_test

        time_p = modelPKN_tp(length_test, q_ft3, frac_ht, fl_ht, w_avg_in, cl, spurt_ft)
        # Vi in bbl
        Vi_test = Q * time_p # required Vi

        eff_tmp = (A_f * w_avg_in / 12 / 5.61) / Vi_test
        alpha = exponentAlfa(n, eff_tmp, "PKN")

        if abs(Vi - Vi_test) < 0.1: # error is less than 0.1
            break
        else:
            length_test *= (Vi / Vi_test)

    return round(length_test, 0)
    
