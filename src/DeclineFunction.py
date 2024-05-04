import math

PI = 3.1416

def asin(x):
    if x == 1:
        asin = PI / 2
    else:
        asin = math.atan(x / math.sqrt(-x * x + 1))
    return asin

def small_g(tD, n_prime_rheo: float = 0.5, estEfficiency: float = 0.5, model2D: str = 'PKN'):
    alfa = exponentAlfa(n_prime_rheo, estEfficiency, model2D)
    if tD > 0:
        small_g = g_up(tD) - 2 * (1 - alfa) * (g_up(tD) - g_lw(tD))
    else:
        small_g = (4 / 3 - 2 * (1 - alfa) * (4 / 3 - PI / 2))
    return small_g

def g_up(tD):
    g_up = 4 / 3 * (((1.0 + tD) ** 1.5) - tD ** 1.5)

    return g_up

def g_lw(self,tD):
    g_lw = ((1.0 + tD) * asin((1.0 + tD) ** -0.5)) + math.sqrt(tD)

    return g_lw

def F_L(t_fromStartPumping, tc):
    t_ov_tc = t_fromStartPumping / tc
    
    if t_ov_tc > 1:
        f_l = 2 / PI * asin(math.sqrt(1 / t_ov_tc))
    else:
        return "Error."

    return f_l

def logLogSlope(x, y, offset):
    X1 = math.log(x.offset(-offset, 0))

def dydx(X1, x, X2, Y1, y, Y2):
    dx1 = x - X1
    dx2 = X2 - x
    dy1 = y - Y1
    dy2 = Y2 - y

    if X2 != x and x != X1:
        dydx = abs((dx2 * (dy1 / dx1) + dx1 * (dy2 / dx2)) / (dx1 + dx2))
    else:
        return "Error."
    
    return dydx

def G(tD_shutin, n_prime_rheo: float = 0.5, estEfficiency: float = 0.5, model2D: str = 'PKN'):
    tD = tD_shutin
    if tD >= 0:
        g0 = small_g(0, n_prime_rheo, estEfficiency, model2D)
        g_tD = small_g(tD, n_prime_rheo, estEfficiency, model2D)
        G = 4 / PI * (g_tD - g0)
    else:
        return "Error."
    
    return G

def exponentAlfa(n_prime, eff, model2D):
    if model2D == 'PKN':
        alpha_up = (2.0 * n_prime + 2.0) / (2.0 * n_prime + 3.0)
    elif model2D == 'KGD':
        alpha_up = (n_prime + 1.0) / (n_prime + 2.0)
    elif model2D == 'RADIAL':
        alpha_up = (4.0 * n_prime + 4.0) / (3.0 * n_prime + 6.0)

    alpha_lw = 0.5

    return eff * alpha_up + (1.0 - eff) * alpha_lw