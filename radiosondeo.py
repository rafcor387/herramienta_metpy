import metpy.calc as mpcalc
from metpy.units import units
from metpy.interpolate import interpolate_1d, log_interpolate_1d
import metpy.constants as mpconst
import pandas as pd
import numpy as np

#niveles
_P = list(np.arange(620.0, 80.0, -20.0, dtype=float))
P_LEVELS = np.array(_P, dtype=float)*units.hPa   
P_LEVELS_ASC = np.sort(P_LEVELS)

""" si la pressure >=100  debe tomar esa linea y llegar hasta ahi en todos los radiosondeos"""
df = pd.read_csv("Radiosondas-2018/02022018EDT.tsv", sep="\t", skiprows=45)
df.columns = df.columns.str.strip()
df["P"] = df["P"].astype(float)
df = df[df['P'] >= 100].copy()

#principales 
pressure = df['P'].to_numpy() * units.hPa
height = df['Height'].to_numpy() * units.meter
temperature = df['T'].to_numpy() * units.kelvin
dewpoint = df['TD'].to_numpy() * units.kelvin
wind_speed = df['FF'].to_numpy() * units.meter / units.second
wind_dir = df['DD'].to_numpy() * units.degrees

tempC = (df['T'] - 273.15).to_numpy() * units.degC

#secundarias
v = df['v'].to_numpy() *units.meter / units.second
u = df['u'].to_numpy() *units.meter / units.second
rh = df['RH'].to_numpy() * units.percent
mr = df['MR'].to_numpy() * units('g/kg')

#CALCULOS 
#viento 
uc, vc = mpcalc.wind_components(wind_speed, wind_dir) 
# ========= TERMODINÁMICA BÁSICA =========
# Humedad y termodinámica
#rh = mpcalc.relative_humidity_from_dewpoint(temperature, dewpoint)                # 0-1
e  = mpcalc.vapor_pressure(pressure, mr)                                  # presión de vapor
es = mpcalc.saturation_vapor_pressure(temperature)                           # presión de vapor saturado
ws = mpcalc.saturation_mixing_ratio(pressure, temperature)                          # razón de mezcla de saturación
w  = mpcalc.mixing_ratio_from_relative_humidity(pressure, temperature, rh)          # razón de mezcla
Tv = mpcalc.virtual_temperature(temperature, mr, mpconst.nounit.epsilon)                              # T virtual
theta = mpcalc.potential_temperature(pressure, temperature)                         # θ
theta_v = mpcalc.virtual_potential_temperature(pressure, temperature, mr, mpconst.epsilon)            # θv
theta_e = mpcalc.equivalent_potential_temperature(pressure, temperature, dewpoint)        # θe (versión general)
rho = mpcalc.density(pressure, temperature, mr, mpconst.epsilon)                                         # kg/m^3 

# ========= NIVELES CARACTERÍSTICOS Y PARCELAS =========
# Parcela de superficie (SB)
p_sfc, T_sfc, Td_sfc = pressure[0], temperature[0], dewpoint[0]
prof_sb = mpcalc.parcel_profile(pressure, T_sfc, Td_sfc)
lcl_p, lcl_T = mpcalc.lcl(p_sfc, T_sfc, Td_sfc)
lfc_p, lfc_T = mpcalc.lfc(pressure, temperature, dewpoint, prof_sb, Td_sfc)
el_p, el_T  = mpcalc.el(pressure, temperature, dewpoint, prof_sb)

cape_sb, cin_sb = mpcalc.cape_cin(pressure, temperature, dewpoint, prof_sb)

# Parcela mixta mixed_parcel (lowest 100 hPa por defecto; ajusta si quieres)
p_mp, T_mp, Td_mp = mpcalc.mixed_parcel(pressure, temperature, dewpoint, depth=50 * units.hPa)
# Capa_mixta mixed_layer 
#prof_ml = mpcalc.parcel_profile(p_ml, T_ml, Td_ml)
cape_ml, cin_ml = mpcalc.mixed_layer_cape_cin(pressure, temperature, dewpoint, depth=50 * units.hPa)
#mixed_layer

""" # ========= ÍNDICES DE ESTABILIDAD CLÁSICOS =========
k_index = mpcalc.k_index(pressure, temperature, dewpoint)
showalter = mpcalc.showalter_index(pressure, temperature, dewpoint)
lifted = mpcalc.lifted_index(pressure, temperature, prof_sb)  # LI respecto a 500 hPa
tt_index = mpcalc.total_totals_index(pressure, temperature, dewpoint)    
 """
#cape_cin, surface_based_cape_cin, most_unstable_cape_cin, mixed_layer_cape_cin

#interpolacion
# Interpolamos temperatura y punto de rocío a los nuevos niveles
temp_interp = log_interpolate_1d(P_LEVELS_ASC, pressure, temperature)
dew_interp = log_interpolate_1d(P_LEVELS_ASC, pressure, dewpoint)

print("Datos principales")
print("presion:", pressure)
print("altura:", height)
print("temperatura en kelvin:", temperature)
print("temperatura de rocío en kelvin:", dewpoint)
print("velocidad del viento:", wind_speed)
print("dirección del viento:", wind_dir)
print("temperatura en centigrados:", tempC)

print("\nDatos secundarios")
print("componente v:", v)
print("componente u:", u)
print("humedad relativa:", rh)
print("mezcla de razón:", mr)

print("\nCálculos") 
print("componentes del viento (u)", uc) 
print("componentes del viento (v)", vc) 

print("\nInterpolación:", P_LEVELS_ASC) 
print("temperature:", temp_interp) 
print("dewpoint:", dew_interp) 

print("\nHumedad relativa calculada:", rh) 
print("presión de vapor", e)
print("mezcla de razón ws:", ws)
print("mezcla de razón w:", w)

print("\n========= VARIABLES CALCULADAS =========")
print("Presión de vapor (e):", e)
print("Presión de vapor saturado (es):", es)
print("Razón de mezcla de saturación (ws):", ws)
print("Razón de mezcla (w):", w)
print("Temperatura virtual (Tv):", Tv)
print("Temperatura potencial (θ):", theta)
print("Temperatura potencial virtual (θv):", theta_v)
print("Temperatura potencial equivalente (θe):", theta_e)
print("Densidad del aire (ρ):", rho)

print("\n========= NIVELES Y PARCELAS =========")
print("Presión superficie (Psfc):", p_sfc)
print("Temperatura superficie (Tsfc):", T_sfc)
print("Temperatura de rocío superficie (Tdsfc):", Td_sfc)
print("Perfil de la parcela de SB (prof_sb):", prof_sb)
print("Nivel de condensación por elevación (LCL) presión:", lcl_p)
print("Nivel de condensación por elevación (LCL) temperatura:", lcl_T)
print("Nivel de libre convección (LFC) presión:", lfc_p)
print("Nivel de equilibrio (EL) presión:", el_p)
print("CAPE superficie (SB):", cape_sb)
print("CIN superficie (SB):", cin_sb)
print("mixed_parcel presion:", p_mp," temperatura:",T_mp, "dewpoint:",Td_mp )
print("cape y cin de mixed_layer:", cape_ml, cin_ml)


def print_var(name, var):
    """Imprime el nombre, unidades y primeros valores de una variable con unidades."""
    try:
        vals = np.asarray(var)
        units_str = str(var.units)
        if vals.ndim == 0:
            print(f"{name:>15}: {vals:.3f} {units_str}")
        else:
            sample = ", ".join([f"{v:.3f}" for v in vals[:5]])
            print(f"{name:>15}: [{sample}, ...] {units_str} (n={len(vals)})")
    except Exception as e:
        print(f"{name:>15}: Error al mostrar ({e})")

