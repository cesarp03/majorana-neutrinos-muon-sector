# Actualizacion de la cota con los parametros de MINOS+ (Adamson et al. 2020)
# y NuFIT 5.2, para las jerarquias normal (NH) e invertida (IH), con analisis
# de sensibilidad a 1 y 2 sigma.
# Autor: Cesar Elias Pena Aparicio
#
# Uso:     python scripts/02_actualizacion_minos_plus.py
# Salida:  tres figuras en figuras/ y resumen por consola.

from pathlib import Path

import numpy as np
import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches

FIGURAS = Path(__file__).resolve().parents[1] / "figuras"
FIGURAS.mkdir(exist_ok=True)

# Configuración fenomenológica (MINOS+ / NuFIT 5.2)
E_nu_GeV = 2.0
L_km = 735.0
lhs_limit = 0.001 * (E_nu_GeV * 1e9)**2

# Parámetros fijos globales
dm21_sq = 7.53e-5
theta12 = np.arcsin(np.sqrt(0.307))
theta13 = np.arcsin(np.sqrt(0.022))

# Valores centrales y errores (adamson2020)
# NH: Normal Hierarchy, IH: Inverted Hierarchy
params_central = {
    'NH': {'dm32_sq': 2.40e-3,  'theta23': np.arcsin(np.sqrt(0.43))},
    'IH': {'dm32_sq': -2.40e-3, 'theta23': np.arcsin(np.sqrt(0.43))}
}

# Motor de cálculo vectorizado
def calcular_cota_vectorizada(a2, a3, delta_cp, dm32_sq, theta23, jerarquia):
    c12, s12 = np.cos(theta12), np.sin(theta12)
    c13, s13 = np.cos(theta13), np.sin(theta13)
    c23, s23 = np.cos(theta23), np.sin(theta23)

    # 1. Matriz PMNS con delta_cp
    U_std = np.array([
        [c12*c13, s12*c13, s13*np.exp(-1j*delta_cp)],
        [-s12*c23 - c12*s23*s13*np.exp(1j*delta_cp), c12*c23 - s12*s23*s13*np.exp(1j*delta_cp), s23*c13],
        [ s12*s23 - c12*c23*s13*np.exp(1j*delta_cp), -c12*s23 - s12*c23*s13*np.exp(1j*delta_cp), c23*c13]
    ])
    P = np.diag([1.0, np.exp(1j*a2), np.exp(1j*a3)])
    U = U_std @ P
    U_mu = U[1, :]

    # 2. Vectorización del rango de masas (Rango extendido: 10^-5 eV a 10^8 eV)
    m_L = np.logspace(-5, 8, 1500) 
    
    if jerarquia == 'NH':
        m1 = m_L
        m2 = np.sqrt(m1**2 + dm21_sq)
        m3 = np.sqrt(m2**2 + dm32_sq)
        DM2 = np.array([[0, -dm21_sq, -(dm21_sq + dm32_sq)],
                        [0, 0, -dm32_sq],
                        [0, 0, 0]])
    else: # IH
        m3 = m_L
        m2 = np.sqrt(m3**2 + abs(dm32_sq))
        m1 = np.sqrt(m2**2 - dm21_sq)
        DM2 = np.array([[0, -dm21_sq, abs(dm32_sq) - dm21_sq],
                        [0, 0, abs(dm32_sq)],
                        [0, 0, 0]])
        
    masses = np.array([m1, m2, m3])
    m_eff_eV = np.abs(U_mu[0]**2 * m1 + U_mu[1]**2 * m2 + U_mu[2]**2 * m3)

    # 3. Cálculo de términos oscilatorios
    osc_term = np.zeros_like(m_L)
    for i in range(3):
        for j in range(i+1, 3):
            U2i_U2jstar = U_mu[i]**2 * np.conj(U_mu[j])**2
            
            # Argumentos ajustados exactamente a la metodología (Factor 1.267)
            arg_base = DM2[i, j] * L_km / E_nu_GeV
            sin2_arg = np.sin(1.267 * arg_base)**2
            sin_arg  = np.sin(2.534 * arg_base)
            
            osc_term += 4 * np.real(U2i_U2jstar) * masses[i] * masses[j] * sin2_arg
            osc_term += 2 * np.imag(U2i_U2jstar) * masses[i] * masses[j] * sin_arg

    rhs = m_eff_eV**2 - osc_term

    # 4. Búsqueda de saturación (rhs > lhs)
    rupturas = rhs > lhs_limit
    if np.any(rupturas):
        idx_ruptura = np.argmax(rupturas)
        return m_eff_eV[idx_ruptura] / 1e6 # Devolver cota en MeV
    else:
        return np.nan # Punto ciego (divergencia)

# Rutina principal: barrido y sensibilidad
def ejecutar_barrido(dm32, theta23, jerarquia):
    res_alpha = 40 # Resolución de fases de Majorana
    res_delta = 10 # Exploración sistemática de Dirac (delta_CP)
    
    alpha_vals = np.linspace(0, 2*np.pi, res_alpha)
    delta_vals = np.linspace(0, 2*np.pi, res_delta)
    
    mapa_cotas = np.zeros((res_alpha, res_alpha))
    
    for i, a2 in enumerate(alpha_vals):
        for j, a3 in enumerate(alpha_vals):
            cota_max_delta = 0.0 # Buscamos la más conservadora respecto a delta_CP
            es_punto_ciego = True
            
            for d_cp in delta_vals:
                cota = calcular_cota_vectorizada(a2, a3, d_cp, dm32, theta23, jerarquia)
                if not np.isnan(cota):
                    es_punto_ciego = False
                    if cota > cota_max_delta:
                        cota_max_delta = cota
            
            mapa_cotas[i, j] = np.nan if es_punto_ciego else cota_max_delta
            
    return alpha_vals, mapa_cotas

# Ejecución paso a paso de la metodología

print("=== INICIANDO FASE 2: ACTUALIZACIÓN MINOS+ Y SENSIBILIDAD ===")
resultados_globales = {}

# Definición de las perturbaciones a +/- 1 sigma y +/- 2 sigma
# (Basado en duplicar los errores de Adamson et al. 2020)
variaciones = {
    'NH': {
        '+1sig': {'dm32_sq': 2.40e-3 + 0.08e-3, 'theta23': np.arcsin(np.sqrt(0.43 + 0.20))},
        '-1sig': {'dm32_sq': 2.40e-3 - 0.09e-3, 'theta23': np.arcsin(np.sqrt(0.43 - 0.04))},
        # Agregamos las variaciones a 2 sigma (doble de la barra de error)
        '+2sig': {'dm32_sq': 2.40e-3 + 0.16e-3, 'theta23': np.arcsin(np.sqrt(0.43 + 0.40))},
        '-2sig': {'dm32_sq': 2.40e-3 - 0.18e-3, 'theta23': np.arcsin(np.sqrt(0.43 - 0.08))}
    },
    'IH': {
        '+1sig': {'dm32_sq': -2.40e-3 - 0.07e-3, 'theta23': np.arcsin(np.sqrt(0.43 + 0.20))},
        '-1sig': {'dm32_sq': -2.40e-3 + 0.08e-3, 'theta23': np.arcsin(np.sqrt(0.43 - 0.04))},
        # Agregamos las variaciones a 2 sigma para Jerarquía Invertida
        '+2sig': {'dm32_sq': -2.40e-3 - 0.14e-3, 'theta23': np.arcsin(np.sqrt(0.43 + 0.40))},
        '-2sig': {'dm32_sq': -2.40e-3 + 0.16e-3, 'theta23': np.arcsin(np.sqrt(0.43 - 0.08))}
    }
}

for jerarquia in ['NH', 'IH']:
    print(f"\n" + "="*50)
    print(f" PROCESANDO JERARQUÍA: {jerarquia}")
    print("="*50)
    
    dm32_c = params_central[jerarquia]['dm32_sq']
    t23_c = params_central[jerarquia]['theta23']
    
    # 1. Caso Central
    A_vals, mapa_central = ejecutar_barrido(dm32_c, t23_c, jerarquia)
    validas_central = mapa_central[~np.isnan(mapa_central)]
    
    print("\n[ VALORES CENTRALES ]")
    print(f"  -> Cota Estricta:     {np.min(validas_central):.2f} MeV")
    print(f"  -> Cota Mediana:      {np.median(validas_central):.2f} MeV")
    print(f"  -> Cota Conservadora: {np.max(validas_central):.2f} MeV")
    
    resultados_globales[jerarquia] = mapa_central
    
    # Análisis de sensibilidad a 1 y 2 sigma
    for sig_label in ['+1sig', '-1sig', '+2sig', '-2sig']:
        print(f"\n[ SENSIBILIDAD {sig_label} ]")
        dm32_v = variaciones[jerarquia][sig_label]['dm32_sq']
        t23_v = variaciones[jerarquia][sig_label]['theta23']
        
        _, mapa_var = ejecutar_barrido(dm32_v, t23_v, jerarquia)
        validas_var = mapa_var[~np.isnan(mapa_var)]
        
        if len(validas_var) > 0:
            print(f"  -> Cota Estricta:     {np.min(validas_var):.2f} MeV")
            print(f"  -> Cota Mediana:      {np.median(validas_var):.2f} MeV")
            print(f"  -> Cota Conservadora: {np.max(validas_var):.2f} MeV")
        else:
            print("  -> Todo el espacio divergió (Punto ciego global).")

# Mapas de calor para las jerarquías normal e invertida
fig, axes = plt.subplots(1, 2, figsize=(14, 6))
cmap = plt.cm.magma.copy()
cmap.set_bad(color='gray')

A2, A3 = np.meshgrid(A_vals, A_vals)

for ax, h in zip(axes, ['NH', 'IH']):
    im = ax.pcolormesh(A2, A3, resultados_globales[h].T, cmap=cmap, shading='auto')
    ax.set_title(f'Espacio de Fases - MINOS+ 2020 ({h})')
    ax.set_xlabel(r'Fase $\alpha_2$ (rad)')
    ax.set_ylabel(r'Fase $\alpha_3$ (rad)')
    fig.colorbar(im, ax=ax, label=r'Cota $m_{eff}$ (MeV)')

gray_patch = mpatches.Patch(color='gray', label='Puntos ciegos')
fig.legend(handles=[gray_patch], loc='upper center')

plt.tight_layout(rect=[0, 0, 1, 0.95])
fig.savefig(FIGURAS / "02_mapas_minos_plus_NH_IH.png", dpi=150)

# Calcular la diferencia absoluta entre mapas
# NOTA: np.nan - valor sigue siendo np.nan, lo cual respeta los puntos ciegos
mapa_diferencia = resultados_globales['NH'] - resultados_globales['IH']

plt.figure(figsize=(7, 6))

# Encontramos el valor máximo absoluto para centrar el mapa en 0
max_val = np.nanmax(np.abs(mapa_diferencia))

# Graficamos con un mapa divergente centrado en 0
im = plt.pcolormesh(A2, A3, mapa_diferencia.T, 
                    cmap='coolwarm', 
                    vmin=-max_val, vmax=max_val, 
                    shading='auto')

plt.title('Diferencia Neta de Cotas: NH - IH')
plt.xlabel(r'Fase $\alpha_2$ (rad)')
plt.ylabel(r'Fase $\alpha_3$ (rad)')
cbar = plt.colorbar(im)
cbar.set_label(r'$\Delta m_{eff}$ (MeV)  [Rojo: NH > IH  |  Azul: IH > NH]')
plt.grid(True, linestyle='--', alpha=0.5)
plt.tight_layout()
plt.savefig(FIGURAS / "02_diferencia_NH_IH.png", dpi=150)


fig, ax = plt.subplots(figsize=(8, 6))

# Fondo: Mapa de calor de NH
im = ax.pcolormesh(A2, A3, resultados_globales['NH'].T, cmap='magma', alpha=0.8)

# Cinco niveles de masa equiespaciados entre el mínimo y el máximo de NH
niveles = np.linspace(np.nanmin(resultados_globales['NH']), np.nanmax(resultados_globales['NH']), 5)

# Contornos de NH (Líneas continuas blancas)
cs_nh = ax.contour(A2, A3, resultados_globales['NH'].T, levels=niveles, colors='white', linewidths=1.5)
ax.clabel(cs_nh, inline=True, fmt='%.1f', fontsize=9)

# Contornos de IH (Líneas discontinuas cian)
cs_ih = ax.contour(A2, A3, resultados_globales['IH'].T, levels=niveles, colors='cyan', linestyles='dashed', linewidths=1.5)

# Leyenda de los contornos
h1, _ = cs_nh.legend_elements()
h2, _ = cs_ih.legend_elements()
ax.legend([h1[0], h2[0]], ['Contornos NH', 'Contornos IH'], loc='lower left')

plt.title('Desplazamiento Geométrico de Cotas (NH vs IH)')
plt.tight_layout()
plt.savefig(FIGURAS / "02_contornos_NH_IH.png", dpi=150)
