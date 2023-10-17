import numpy as np
from numpy.typing import NDArray
import matplotlib.pyplot as plt
import customtkinter as ctkinter

class hydraulicFracture:
  def __init__(self, model: str = 'PKN'):
    self.model = model

    self.w_avg_ft = self.w_avg_in / 12
    self.spurt_ft = self.spurt * 0.001336806 # convert from gal/100ft2 to ft
    self.q_ft3 = self.q * 5.614583

  def frac_length(self):
    """
    Find fracture length based on PKN model
    """
    self.plain_strain = self.young / (1 - self.poisson**2)
    k_apos = self.app_visc / (47880 * (100**(self.n_shear-1)))

    length_test = 500 # in ft
    eff_tmp = 0.5

    for _ in range(1,50):
      w_max = 3 * ((k_apos * length_test) * (self.frac_height / self.plain_strain) * (self.q_ft3 / (self.frac_height * 60))**self.n) ** (1 / (2 * self.n+2))
      beta_p = (self.n + 2) / (self.n + 3)
      self.w_avg_in = 12 * self.PI / 4 * beta_p * w_max
      A_f = 2 * self.frac_height * length_test

      time_p = self.modelPKN_tp(length_test, self.q_ft3, self.frac_height, self.fl_height, self.w_avg_in, self.cl, self.spurt_ft)
      Vi_test = self.q * time_p

      eff_tmp = (A_f * self.w_avg_in / 12 / 5.61) / Vi_test
      self.alpha = self.exp_alpha(self.n, eff_tmp, self.model)

      if np.abs(self.Vi - Vi_test) < 0.1:
        break
      else:
        length_test *= (self.Vi / Vi_test)

    return np.round(length_test, 0)
  
  def modelPKN_tp(self, frac_length, q, frac_height, fl_height, w_avg_in, cl, spurt_ft):
    """
    Calculates the pumping time (tp)
    """
    r_p = fl_height / frac_height
    A_f = 2 * frac_height * frac_length

    go = 1.5
    self.vol_lost = 2 * r_p * A_f * go * cl

    a = q
    b = -3 * A_f * r_p * cl
    c = -A_f * (w_avg_in + 2 * r_p * spurt_ft)

    sqr_t = (-b + np.sqrt(b**2 - 4 * a * c)) / (2 * a)
    
    return sqr_t**2
  
  def design_length(self):
    """
    This one is p weird lmao
    """
    w_max = 3 * 12 * ((self.k_apos * self.frac_length) * (self.frac_height / self.plain_strain) *(self.q_ft3 / (self.frac_height * 60))**self.n) ** (1 / (2 * self.n+2))
    net_press = (self.young * w_max) / (2 * self.frac_height)
    beta_p = (self.n + 2) / (self.n + 3)
    w_avg = np.pi * beta_p * w_max
    
    return net_press, w_avg
  
  def pump_time(self):
    return
  
  def exp_alpha(self, n_prime, eff, model):
    """
    Calculates exponential alpha... idk what the fuck this is for (yet)

    Args
    ----------------------
    n_prime: float
      ???
    eff: float
      Efficiency?
    model: string
      Selection of models, PKN, KGD, or RADIAL 2D model.

    Returns
    ----------------------
    The exponential alpha ^^
    """
    if model == 'PKN':
      self.alpha_up = (2.0 * n_prime + 2.0) / (2.0 * n_prime + 3.0)
    elif model == 'KGD':
      self.alpha_up = (n_prime + 1.0) / (n_prime + 2.0)
    elif model == 'RADIAL':
      self.alpha_up = (4.0 * n_prime + 4.0) / (3.0 * n_prime + 6.0)

    alpha_lw = 0.5

    return eff * self.alpha_up + (1.0 - eff) * alpha_lw

  
class appInput:
  def __init__(self):
    self.geometry("600x500")
    self.title("Hydraulic Fracturing 101")

    self.button = ctkinter.CTkButton(self, text=f'Welcome to hydraulic fracturing app. You have chosen the {hydraulicFracture.model} model. Please input the data by clicking on the button below.', command=self.button_click)
    self.button.grid(row=0, column=0, padx=20, pady=10)

  def button_click(self):
    dialog = ctkinter.CTkInputDialog(text='')