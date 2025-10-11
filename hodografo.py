import matplotlib.pyplot as plt
import numpy as np
from metpy.plots import Hodograph
from metpy.units import units
import pandas as pd

df = pd.read_csv("./Radiosondas-2018/02022018EDT.tsv", sep="\t", skiprows=44)
df.columns = df.columns.str.strip()

pressure = df['P'].values * units.hPa

v_wind = df['v'].values * units.meter / units.second
u_wind = df['u'].values * units.meter / units.second


min_len = min(len(pressure), len(u_wind), len(v_wind))
pressure = pressure[:min_len]
u_wind = u_wind[:min_len]
v_wind = v_wind[:min_len]


fig = plt.figure(figsize=(6, 6)) 
ax_hodo = fig.add_subplot(1, 1, 1) 

hodo = Hodograph(ax_hodo, component_range=60.) 


hodo.add_grid(increment=5) 
hodo.add_grid(increment=10, ls='-') 

hodo.plot_colormapped(u_wind, v_wind, pressure, cmap='viridis')

ax_hodo.set_xlabel('U-wind (m/s)')
ax_hodo.set_ylabel('V-wind (m/s)')
ax_hodo.set_title('Hodograph Example', fontsize=14)

plt.show()