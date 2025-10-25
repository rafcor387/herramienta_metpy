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

df = pd.read_csv("Radiosondas-2018/20181128EDT.tsv", sep="\t", skiprows=45)
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
a= mpcalc.precipitable_water
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
cape_sb, cin_sb = mpcalc.surface_based_cape_cin(pressure, temperature, dewpoint)

# Parcela mixta (mp) (lowest 100 hPa por defecto; ajusta si quieres)
p_mp, T_mp, Td_mp = mpcalc.mixed_parcel(pressure, temperature, dewpoint, depth=50 * units.hPa)
prof_mp = mpcalc.parcel_profile(pressure,T_mp,Td_mp)
cape_ml, cin_ml = mpcalc.mixed_layer_cape_cin(pressure, temperature, dewpoint, depth=50 * units.hPa)

#Parcel mas inestable (mu)
cape_mu, cin_mu = mpcalc.most_unstable_cape_cin(pressure,temperature, dewpoint)

lcl_p, lcl_T = mpcalc.lcl(p_sfc, T_sfc, Td_sfc)
lfc_p, lfc_T = mpcalc.lfc(pressure, temperature, dewpoint)
el_p, el_T  = mpcalc.el(pressure, temperature, dewpoint)
ccl_1, ccl_2, ccl_3 = mpcalc.ccl(pressure, temperature, dewpoint)

# ========= ÍNDICES DE ESTABILIDAD CLÁSICOS =========
lifted = mpcalc.lifted_index(pressure, temperature, prof_sb)  # LI respecto a 500 hPa 

# ========= LAPSE RATE =========
gama_most = mpcalc.moist_lapse(pressure, temperature[0]).to('degC')
dT_dz = (-mpcalc.first_derivative(temperature, x=height))
ELR   = (-dT_dz).to('degC/km')  

# 1) Construir altura AGL (0 en superficie del sondeo)
agl = (height - height[0]).to('meter')
# 2) Extraer exactamente la capa 0–3 km AGL (interpola bordes si no caen en un nivel exacto)
p_03, T_03, agl_03 = mpcalc.get_layer(pressure, temperature, agl, height=agl,
                                      bottom=0*units.meter, depth=3000*units.meter)
# 3) Lapse rate medio en la capa (ambiental):
Gamma_03 = ((T_03[0] - T_03[-1]) / (agl_03[-1] - agl_03[0])).to('degC/km')
#print(Gamma_03)  # → °C/km

p_36, T_36, agl_36 = mpcalc.get_layer(pressure, temperature, agl, height=agl,
                                      bottom=3000*units.meter, depth=3000*units.meter)
Gamma_36 = ((T_36[0] - T_36[-1]) / (agl_36[-1] - agl_36[0])).to('degC/km')


sigma = mpcalc.static_stability(pressure,temperature)
z = mpcalc.height_to_geopotential(height)

#most_estatic_energy = mpcalc.moist_static_energy(height,temperature, humedadEspecifica)   dry_static_energy


#interpolacion
# Interpolamos temperatura y punto de rocío a los nuevos niveles
""" temp_interp = log_interpolate_1d(P_LEVELS_ASC, pressure, temperature)
dew_interp = log_interpolate_1d(P_LEVELS_ASC, pressure, dewpoint) """

print("Datos principales")
print("presion:", pressure[:10])
print("altura:", height[:10])
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

""" print("\nInterpolación:", P_LEVELS_ASC) 
print("temperature:", temp_interp) 
print("dewpoint:", dew_interp)  """

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
print("Nivel de condensación por elevación (LCL) presión:", lcl_p, lcl_T)
print("Nivel de libre convección (LFC) presión:", lfc_p, lfc_T)
print("Nivel de equilibrio (EL) presión:", el_p, el_T)
print("Nivel de condensacion convectiva (CCL):", ccl_1, ccl_2, ccl_3)

print("Presión_sb, temperatura_sb, dewpoint_sb:", p_sfc, T_sfc, Td_sfc)
print("Perfil de la parcela de SB (prof_sb):", prof_sb)
print("parcel por mp:", prof_mp)

print("Presion_mp, temperatura_mp, dewpoint_mp:", p_mp, T_mp,Td_mp )
print("CAPE y CIN superficie (SB):", cape_sb, cin_sb)
print("CAPE y CIN mixed_layer:", cape_ml, cin_ml)
print("CAPE y CIN most_unstable:", cape_mu, cin_mu)

print("\n========= INDICES DE ESTABILIDAD CLASICOS =========")
print("lifted index value:", lifted)

print("\n========= LAPSE RATE (GAMMA) =========")
print("lapse rate moist:", gama_most[:5])
print("lapse ambiental:" , ELR[:5])
print("lapse rate de 0-3 km:", Gamma_03)
print("lapse rate de 3-6 km:", Gamma_36)

print("\n========= ESTABILIDAD ESTATICA (SIGMA) =========")
print("estabilidad estatica:", sigma[:10])

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

