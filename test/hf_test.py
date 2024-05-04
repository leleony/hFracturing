import math
import pandas as pd
import json
import os
import logging

inputDir = './hFracturing/test/data/'

logging.info("Importing data...")

try:
    csvData = pd.read_csv(inputDir+'timestep_data.csv').values
    with open('static_data.json') as files:
        jsonData = json.loads(files)
except:
    logging.error("There's no available data, or it doesn't match the required extension (.csv/.json)!")

# getting available data
## very minimal input
E = jsonData['young']
v = jsonData['poisson']
frac_height = jsonData['frac_height']
fluid_loss_height = jsonData['fl_height']
q_inj = jsonData['q_inj']

plain_strain_mod = E / (1 - v**2)

## frac fluid rheology
n_shear_frac = jsonData["n'"]
app_visc = jsonData['app_visc']

# tip screen out input
add_tp = jsonData["additional_pump_time"]