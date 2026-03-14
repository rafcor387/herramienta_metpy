import pandas as pd
import numpy as np
import metpy.calc as mpcalc
from metpy.units import units
import metpy.constants as mpconsts
from scipy.integrate import trapezoid

def cargar_y_limpiar_datos(ruta_archivo):
    """
    Carga el radiosondeo y soluciona el problema de presiones duplicadas
    y perfiles no ordenados.    
    """
    # Usamos skiprows=47 basándonos en la estructura de tus archivos .tsv
    df = pd.read_csv(ruta_archivo, sep="\t", skiprows=45)
    df.columns = df.columns.str.strip()

    #Convertir a numérico y eliminar nulos
    df['P'] = pd.to_numeric(df['P'], errors='coerce')
    df = df.dropna(subset=['P', 'T', 'TD'])
    
    #Eliminar filas con presión duplicada (nos quedamos solo con la primera vez que aparece)
    df = df.drop_duplicates(subset=['P'], keep='first')
    
    # Asegurar que el perfil sea estrictamente decreciente
    df = df.sort_values(by='P', ascending=False)

    # Extraemos variables con sus unidades físicas de MetPy
    p = df['P'].values * units.hPa
    t = df['T'].values * units.kelvin
    td = df['TD'].values * units.kelvin
    h = df['Height'].values * units.meter

    return p, t, td, h

def calcular_termodinamica_parcela(p, t, td):
    """
    Calcula los parámetros principales de la termodinámica de la parcela.
    """
    # 1. Parámetros de la Parcela de Superficie (SB)
    # -----------------------------------------------
    # Nivel de Condensación por Ascenso (LCL)
    lcl_p, lcl_t = mpcalc.lcl(p[0], t[0], td[0])

    # Perfil de la parcela (la curva que seguiría el aire al subir)
    perfil_parcela = mpcalc.parcel_profile(p, t[0], td[0]).to('degK')

    # CAPE y CIN basados en superficie
    cape_sb, cin_sb = mpcalc.surface_based_cape_cin(p, t, td)

    # 2. Parámetros de la Parcela de Capa de Mezcla (ML)
    # -------------------------------------------------
    # En La Paz (aprox 620 hPa), una capa de 50 hPa es más realista
    cape_ml, cin_ml = mpcalc.mixed_layer_cape_cin(p, t, td, depth=50 * units.hPa)

    # 3. Niveles Característicos (LFC y EL)
    # -------------------------------------
    lfc_p, lfc_t = mpcalc.lfc(p, t, td)
    el_p, el_t = mpcalc.el(p, t, td)

    # 4. Contenido de Agua
    # --------------------
    agua_precipitable = mpcalc.precipitable_water(p, td)

    # Estructura de salida
    resultados = {
        "LCL": {"presion": lcl_p, "temperatura": lcl_t},
        "LFC": {"presion": lfc_p, "temperatura": lfc_t},
        "EL": {"presion": el_p, "temperatura": el_t},
        "CAPE_SB": cape_sb,
        "CIN_SB": cin_sb,
        "CAPE_ML": cape_ml,
        "CIN_ML": cin_ml,
        "PW": agua_precipitable,
        "perfil_parcela": perfil_parcela
    }

    return resultados

def analisis_termodinamico_total(p, t, td, h):
    """
    Ejecuta todas las funciones de termodinámica de parcela disponibles en MetPy
    para un perfil completo.
    """
    # --- A. DEFINICIÓN DE PARCELAS ---
    # 1. Parcela de Superficie (SB): Tal cual está el suelo (626 hPa)
    p_sb, t_sb, td_sb = p[0], t[0], td[0]

    # 2. Parcela más Inestable (MU): Busca el aire con más energía en los primeros 300hPa
    p_mu, t_mu, td_mu, indice_mu = mpcalc.most_unstable_parcel(p, t, td, depth=300 * units.hPa)

    # 3. Parcela de Capa de Mezcla (ML): Promedia el aire de los primeros 50hPa (Ideal para La Paz)
    p_ml, t_ml, td_ml = mpcalc.mixed_parcel(p, t, td, depth=50 * units.hPa)

    # --- B. CÁLCULO DE PERFILES (La trayectoria de la burbuja) ---
    prof_sb = mpcalc.parcel_profile(p, t_sb, td_sb).to('degK')
    prof_mu = mpcalc.parcel_profile(p, t_mu, td_mu).to('degK')
    prof_ml = mpcalc.parcel_profile(p, t_ml, td_ml).to('degK')

    # --- C. ENERGÍA (CAPE y CIN) PARA CADA PARCELA ---
    # Aquí calculamos los 3 tipos de CAPE y CIN que existen
    cape_sb, cin_sb = mpcalc.surface_based_cape_cin(p, t, td)
    cape_mu, cin_mu = mpcalc.most_unstable_cape_cin(p, t, td, depth=300 * units.hPa)
    cape_ml, cin_ml = mpcalc.mixed_layer_cape_cin(p, t, td, depth=50 * units.hPa)

    # --- D. NIVELES CRÍTICOS (LCL, LFC, EL, CCL) ---
    # LCL: Nivel de Condensación (Donde empieza la nube si sube por fuerza externa)
    lcl_p, lcl_t = mpcalc.lcl(p_sb, t_sb, td_sb)

    # LFC: Nivel de Libre Convección (Donde la nube flota sola por ser más cálida)
    lfc_p, lfc_t = mpcalc.lfc(p, t, td)

    # EL: Nivel de Equilibrio (Donde la nube deja de subir, el "techo")
    el_p, el_t = mpcalc.el(p, t, td)

    # CCL y Temperatura Convectiva (Tc): ¡Vital!
    # El CCL es donde se forman nubes solo por el calor del sol (sin ayuda de frentes)
    # Tc es la temperatura que debe alcanzar el suelo en La Paz para que empiecen las nubes.
    ccl_p, ccl_t, temp_convectiva = mpcalc.ccl(p, t, td)

    # --- E. PARÁMETROS DE DESCENSO Y HUMEDAD ---
    # DCAPE (Downdraft CAPE): Energía de las corrientes descendentes (Riesgo de micro-ráfagas)
    dcape, d_t = mpcalc.downdraft_cape(p, t, td)

    # Agua Precipitable Total
    pw = mpcalc.precipitable_water(p, td)

    return {
        "SB": {"CAPE": cape_sb, "CIN": cin_sb, "LCL": lcl_p},
        "MU": {"CAPE": cape_mu, "CIN": cin_mu, "LFC": lfc_p, "P_Origen": p_mu},
        "ML": {"CAPE": cape_ml, "CIN": cin_ml, "EL": el_p},
        "CONVECCION_SOLAR": {"CCL": ccl_p, "Temp_Convectiva": temp_convectiva.to('degC')},
        "PELIGROS": {"DCAPE": dcape, "PW": pw},
        "perfiles": {"prof_sb": prof_sb, "prof_mu": prof_mu, "prof_ml": prof_ml}
    }
    
def calcular_termodinamica_segura(p, t, td):
    """
    Calcula termodinámica de parcela adaptándose al límite del sensor.
    """
    # 1. Verificar qué tan alto llegó el globo
    presion_minima_archivo = np.min(p) # El punto más alto alcanzado
    presion_superficie = p[0]
    profundidad_total_disponible = presion_superficie - presion_minima_archivo

    # 2. Ajustar los límites de búsqueda (Depths)
    # Si el globo no llegó muy alto, reducimos la búsqueda para no causar error
    profundidad_mu = min(300 * units.hPa, profundidad_total_disponible * 0.9)
    profundidad_ml = min(50 * units.hPa, profundidad_total_disponible * 0.2)

    # 3. Diccionario para guardar resultados (vacío al inicio)
    res = {}

    # 4. Cálculos con bloques 'try-except' para mayor seguridad
    try:
        # Parcela más inestable (MU) ajustada
        p_mu, t_mu, td_mu, _ = mpcalc.most_unstable_parcel(p, t, td, depth=profundidad_mu)
        cape_mu, cin_mu = mpcalc.most_unstable_cape_cin(p, t, td, depth=profundidad_mu)
        res['MU_CAPE'] = cape_mu
    except Exception as e:
        res['MU_CAPE'] = "Error: Datos insuficientes para MU"

    try:
        # Nivel de Condensación Convectiva (CCL)
        # Esta función a veces falla si el aire es extremadamente seco
        ccl_p, ccl_t, t_conv = mpcalc.ccl(p, t, td)
        res['CCL_P'] = ccl_p
        res['T_CONVECTIVA'] = t_conv.to('degC')
    except Exception as e:
        res['CCL_P'] = "Error: Aire muy seco para CCL"

    # Agua precipitable (esta casi nunca falla si hay datos)
    res['PW'] = mpcalc.precipitable_water(p, td)

    return res

def calcular_dcape_lapaz(pressure, temperature, dewpoint):
    """
    Calcula el DCAPE adaptado para la altitud de La Paz (superficie < 700 hPa).
    """
    try:
        # 1. Definir la capa de origen para el descenso
        # En lugar de 700 hPa, usamos la superficie real de La Paz
        p_base = pressure[0]
        p_tope = 500 * units.hPa

        # Extraemos la capa entre la superficie (~626) y 500 hPa
        p_layer, t_layer, td_layer = mpcalc.get_layer(
            pressure, temperature, dewpoint,
            bottom=p_base,
            depth=p_base - p_tope,
            interpolate=True
        )

        # 2. Hallar la parcela con la energía mínima (Theta-e mínima)
        theta_e = mpcalc.equivalent_potential_temperature(p_layer, t_layer, td_layer)
        idx_min = np.argmin(theta_e)
        p_inicio = p_layer[idx_min]
        t_inicio = t_layer[idx_min]
        td_inicio = td_layer[idx_min]

        # 3. Calcular la temperatura de bulbo húmedo en el inicio del descenso
        t_bulbo_humedo = mpcalc.wet_bulb_temperature(p_inicio, t_inicio, td_inicio)

        # 4. Trayectoria de descenso adiabático húmedo hasta la superficie
        presiones_descenso = pressure[pressure >= p_inicio].to(units.hPa)
        perfil_descenso = mpcalc.moist_lapse(
            presiones_descenso, t_bulbo_humedo, reference_pressure=p_inicio
        )

        # 5. Temperaturas virtuales para el cálculo de flotabilidad negativa
        # Parcela (descendiendo saturada)
        t_v_parcela = mpcalc.virtual_temperature_from_dewpoint(presiones_descenso, perfil_descenso, perfil_descenso)
        # Entorno (aire real medido por el globo)
        t_v_entorno = mpcalc.virtual_temperature_from_dewpoint(
            presiones_descenso,
            temperature[pressure >= p_inicio],
            dewpoint[pressure >= p_inicio]
        )

        # 6. Integración numérica (Método del trapecio como hace MetPy)
        diff = (t_v_entorno - t_v_parcela).to(units.degK).magnitude
        lnp = np.log(presiones_descenso.magnitude)

        dcape_valor = -(mpconsts.Rd * units.Quantity(trapezoid(diff, lnp), 'K')).to(units('J/kg'))

        return dcape_valor

    except Exception as e:
        return f"No se pudo calcular DCAPE: {e}"

# Ejemplo de integración en tu reporte:
# dcape_final = calcular_dcape_lapaz(p, t, td)

def analizar_termodinamica_completa(p, t, td):
    """Extrae absolutamente todos los parámetros de parcela de MetPy."""

    # --- A. LÍMITES DE SEGURIDAD PARA LA PAZ ---
    presion_suelo = p[0]
    presion_tope = np.min(p)
    rango_total = presion_suelo - presion_tope

    # Ajustamos la profundidad de búsqueda para que no se salga del archivo
    prof_mu = min(300 * units.hPa, rango_total * 0.8)
    prof_ml = min(50 * units.hPa, rango_total * 0.2)

    # --- B. CÁLCULOS ---
    resultados = {}

    try:
        # 1. NIVELES CARACTERÍSTICOS
        # LCL: Nivel de Condensación por Ascenso
        lcl_p, lcl_t = mpcalc.lcl(p[0], t[0], td[0])
        resultados['LCL_presion'] = lcl_p

        # LFC: Nivel de Libre Convección
        lfc_p, lfc_t = mpcalc.lfc(p, t, td)
        resultados['LFC_presion'] = lfc_p

        # EL: Nivel de Equilibrio (Techo de la nube)
        el_p, el_t = mpcalc.el(p, t, td)
        resultados['EL_presion'] = el_p

        # CCL: Nivel de Condensación Convectiva y Temperatura Convectiva
        ccl_p, ccl_t, t_conv = mpcalc.ccl(p, t, td)
        resultados['CCL_presion'] = ccl_p
        resultados['Temperatura_Convectiva'] = t_conv.to('degC')

        # 2. ENERGÍA (LAS TRES PARCELAS)
        # CAPE/CIN de Superficie (SB)
        cape_sb, cin_sb = mpcalc.surface_based_cape_cin(p, t, td)
        resultados['SB_CAPE'] = cape_sb
        resultados['SB_CIN'] = cin_sb

        # CAPE/CIN de Capa de Mezcla (ML) - El más real
        cape_ml, cin_ml = mpcalc.mixed_layer_cape_cin(p, t, td, depth=prof_ml)
        resultados['ML_CAPE'] = cape_ml

        # CAPE/CIN Más Inestable (MU) - El peor escenario
        cape_mu, cin_mu = mpcalc.most_unstable_cape_cin(p, t, td, depth=prof_mu)
        resultados['MU_CAPE'] = cape_mu

        # 3. ÍNDICES DE ESTABILIDAD
        # Lifted Index (LI): Estabilidad a 500 hPa
        # Primero necesitamos el perfil de la parcela
        perfil = mpcalc.parcel_profile(p, t[0], td[0])
        li = mpcalc.lifted_index(p, t, perfil)
        resultados['Lifted_Index'] = li[0] # Tomamos el valor a 500hPa

        # 4. CONTENIDO DE AGUA Y DESCENSO
        resultados['Agua_Precipitable'] = mpcalc.precipitable_water(p, td)

        # DCAPE: Energía de las corrientes descendentes
        dcape = calcular_dcape_lapaz(p, t, td)
        resultados['DCAPE'] = dcape

    except Exception as e:
        resultados['Error'] = f"Error en cálculos: {e}"

    return resultados

# --- EJECUCIÓN ---
if __name__ == "__main__":
    archivo = "Radiosondas/01022018EDT.tsv"
    p, t, td, h = cargar_y_limpiar_datos(archivo)
    datos_finales = analizar_termodinamica_completa(p, t, td)

    # Imprimir resultados de forma elegante
    print(f"\n--- REPORTE TERMODINÁMICO: {archivo} ---")
    for clave, valor in datos_finales.items():
        if isinstance(valor, units.Quantity):
            print(f"{clave:25}: {valor:.2f}")
        else:
            print(f"{clave:25}: {valor}")