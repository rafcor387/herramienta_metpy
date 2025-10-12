import pandas as pd
import matplotlib.pyplot as plt 
from metpy.plots import SkewT
import metpy.calc as mpcalc
from metpy.calc import (
    parcel_profile,
    lcl,
    lfc,
    el,
    cape_cin,
    lifted_index,
    wind_components
)
from metpy.units import units

# Leer y limpiar  20181127EDT 20181128EDT.tsv          03052018EDT             29032018EDT       29032018EDT    20181128EDT   09052018EDT
df = pd.read_csv("./Radiosondas-2018/31012018EDT.tsv", sep="\t", skiprows=45)
df.columns = df.columns.str.strip()

# Extraer columnas
pressure = df['P'].values * units.hPa
temperature = df['T'].values * units.kelvin
dewpoint = df['TD'].values * units.kelvin
wind_speed = df['FF'].values * units.meter / units.second
wind_dir = df['DD'].values * units.degrees

# === CÁLCULOS ===
#para la temperatura de la parcela
parcel_prof = parcel_profile(pressure, temperature[0], dewpoint[0])

# CAPE y CIN
cin, cape = cape_cin(pressure, temperature, dewpoint, parcel_prof)
print(f"Skew-T con CAPE={cape:.1f} J/kg, CIN={cin:.1f} J/kg")

#componentes del viento
#u, v = wind_components(wind_speed, wind_dir)

# Niveles importantes
lcl_pressure, lcl_temperature = lcl(pressure[0], temperature[0], dewpoint[0])
lfc_pressure, lfc_temperature = lfc(pressure, temperature, dewpoint)
el_pressure, el_temperature = el(pressure, temperature, dewpoint)

#lifted index
li = lifted_index(pressure, temperature, dewpoint)
print(f"Lifted Index (LI) = {li:.2f}")


# Crear figura
fig = plt.figure(figsize=(9, 9))
skew = SkewT(fig)

# Líneas de referencia
skew.plot_dry_adiabats()
skew.plot_moist_adiabats()
skew.plot_mixing_lines()

# Graficar datos
skew.plot(pressure, temperature, 'red', label="Temp")
skew.plot(pressure, dewpoint, 'green', label="Rocío")
skew.plot(pressure, parcel_prof.to('degK'), 'k', label='Temperatura de la Parcela')
#skew.plot_barbs(pressure[::30], u[::30], v[::30])

# Sombra CAPE y CIN
skew.shade_cin(pressure, temperature, parcel_prof)
skew.shade_cape(pressure, temperature, parcel_prof)

# Marcar niveles importantes
skew.ax.plot(lcl_temperature, lcl_pressure, 'ko', markerfacecolor='cyan', label='LCL')
skew.ax.plot(lfc_temperature, lfc_pressure, 'ko', markerfacecolor='magenta', label='LFC')
skew.ax.plot(el_temperature, el_pressure, 'ko', markerfacecolor='orange', label='EL')

# Ajustes visuales
skew.ax.set_ylim(1000, 100)
skew.ax.set_xlim(-80, 40)
plt.title("Diagrama Skew-T")
plt.legend()
plt.show()

print("valor lfc:", )

# === GUARDAR ETIQUETAS ===""" 
""" # Detectar inversión térmica (opcional)
inversion_mask = temperature[:-1] < temperature[1:]
tiene_inversion = "Sí" if inversion_mask.any() else "No"

# Crear diccionario con las etiquetas
etiquetas = {
    "archivo": "27022018EDT.tsv",  
    "CAPE": round(cape.magnitude, 2),
    "CIN": round(cin.magnitude, 2),
    "Lifted_Index": round(lifted_index(pressure, temperature, dewpoint).magnitude, 2),
    "K_Index": round(mpcalc.k_index(pressure, temperature, dewpoint).magnitude, 2),
    "Total_Totals": round(mpcalc.total_totals_index(pressure, temperature, dewpoint).magnitude, 2),
    "Tiene_inversion": tiene_inversion
}

# Convertir a DataFrame
df_etiquetas = pd.DataFrame([etiquetas])

# Guardar o añadir a archivo CSV
import os

archivo_etiquetas = "etiquetas_radiosondeos.csv"
if not os.path.isfile(archivo_etiquetas):
    df_etiquetas.to_csv(archivo_etiquetas, index=False)
else:
    df_etiquetas.to_csv(archivo_etiquetas, mode='a', header=False, index=False)

print("Etiquetas guardadas en", archivo_etiquetas) """