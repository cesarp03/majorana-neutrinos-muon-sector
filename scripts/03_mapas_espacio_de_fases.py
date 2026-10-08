# Mapas del espacio de fases de Majorana en el regimen de masa pesada (toy model):
# comparacion entre los parametros de 2009 y los de 2020, y entre jerarquias.
# Autor: Cesar Elias Pena Aparicio
#
# Uso:     python scripts/03_mapas_espacio_de_fases.py
# Salida:  tres figuras en figuras/.

from pathlib import Path

import numpy as np
import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches

FIGURAS = Path(__file__).resolve().parents[1] / "figuras"
FIGURAS.mkdir(exist_ok=True)

# Paso 1: configuración del toy model (régimen de masa pesada)
E_NU_GEV = 2.0
L_KM = 735.0
FACTOR_TECNICO = 0.001  
lhs_limit = FACTOR_TECNICO * (E_NU_GEV * 1e9)**2

# Paso 2: motor de cálculo
def calcular_cota_toy_model(a2, a3, delta_cp, dm32_sq, theta23, jerarquia, dm21_sq, theta12, theta13):
    c12, s12 = np.cos(theta12), np.sin(theta12)
    c13, s13 = np.cos(theta13), np.sin(theta13)
    c23, s23 = np.cos(theta23), np.sin(theta23)

    U_std = np.array([
        [c12*c13, s12*c13, s13*np.exp(-1j*delta_cp)],
        [-s12*c23 - c12*s23*s13*np.exp(1j*delta_cp), c12*c23 - s12*s23*s13*np.exp(1j*delta_cp), s23*c13],
        [ s12*s23 - c12*c23*s13*np.exp(1j*delta_cp), -c12*s23 - s12*c23*s13*np.exp(1j*delta_cp), c23*c13]
    ])
    P = np.diag([1.0, np.exp(1j*a2), np.exp(1j*a3)])
    U = U_std @ P
    U_mu = U[1, :]

    m_L = np.logspace(6, 8, 1500) 
    
    if jerarquia == 'NH':
        m1 = m_L
        m2 = np.sqrt(m1**2 + dm21_sq)
        m3 = np.sqrt(m2**2 + abs(dm32_sq))
        DM2 = np.array([[0, -dm21_sq, -(dm21_sq + abs(dm32_sq))], [0, 0, -abs(dm32_sq)], [0, 0, 0]])
    else:
        m3 = m_L
        m2 = np.sqrt(m3**2 + abs(dm32_sq))
        m1 = np.sqrt(m2**2 - dm21_sq)
        DM2 = np.array([[0, -dm21_sq, abs(dm32_sq) - dm21_sq], [0, 0, abs(dm32_sq)], [0, 0, 0]])
        
    masses = np.array([m1, m2, m3])
    m_eff_eV = np.abs(U_mu[0]**2 * m1 + U_mu[1]**2 * m2 + U_mu[2]**2 * m3)

    osc_term = np.zeros_like(m_L)
    for i in range(3):
        for j in range(i+1, 3):
            U2i_U2jstar = U_mu[i]**2 * np.conj(U_mu[j])**2
            arg_base = DM2[i, j] * L_KM / E_NU_GEV
            sin2_arg = np.sin(1.267 * arg_base)**2
            sin_arg  = np.sin(2.534 * arg_base)
            
            osc_term += 4 * np.real(U2i_U2jstar) * masses[i] * masses[j] * sin2_arg
            osc_term += 2 * np.imag(U2i_U2jstar) * masses[i] * masses[j] * sin_arg

    rhs = m_eff_eV**2 - osc_term
    rupturas = rhs > lhs_limit
    
    if np.any(rupturas):
        return m_eff_eV[np.argmax(rupturas)] / 1e6 
    else:
        return np.nan

def ejecutar_barrido_toy(dm32, theta23, jerarquia, dm21, t12, t13, res=40):
    alpha_vals = np.linspace(0, 2*np.pi, res)
    delta_vals = np.linspace(0, 2*np.pi, 10)
    mapa = np.zeros((res, res))
    
    for i, a2 in enumerate(alpha_vals):
        for j, a3 in enumerate(alpha_vals):
            cota_max_delta = -1.0
            es_punto_ciego = True
            for d_cp in delta_vals:
                cota = calcular_cota_toy_model(a2, a3, d_cp, dm32, theta23, jerarquia, dm21, t12, t13)
                if not np.isnan(cota):
                    es_punto_ciego = False
                    if cota > cota_max_delta:
                        cota_max_delta = cota
            mapa[i, j] = np.nan if es_punto_ciego else cota_max_delta
    return alpha_vals, mapa

# Paso 3: mapas con parámetros históricos (2009) y modernos (2020)
dm21_2009, t12_2009, t13_2009 = 7.59e-5, np.arcsin(np.sqrt(0.32)), 0.0  
dm32_2009, t23_2009 = 2.43e-3, np.arcsin(np.sqrt(0.50))

dm21_2020, t12_2020, t13_2020 = 7.53e-5, np.arcsin(np.sqrt(0.307)), np.arcsin(np.sqrt(0.022))
dm32_2020, t23_2020 = 2.40e-3, np.arcsin(np.sqrt(0.43))

print("Calculando mapas del espacio de fases...")
A_vals, mapa_2009_NH = ejecutar_barrido_toy(dm32_2009, t23_2009, 'NH', dm21_2009, t12_2009, t13_2009)
_, mapa_2020_NH = ejecutar_barrido_toy(dm32_2020, t23_2020, 'NH', dm21_2020, t12_2020, t13_2020)
_, mapa_2020_IH = ejecutar_barrido_toy(dm32_2020, t23_2020, 'IH', dm21_2020, t12_2020, t13_2020)

A2, A3 = np.meshgrid(A_vals, A_vals)
cmap_base = plt.cm.magma.copy()
cmap_base.set_bad(color='#d3d3d3') # Gris claro para puntos ciegos
gray_patch = mpatches.Patch(color='#d3d3d3', label='Puntos Ciegos (Interferencia Destructiva Perfecta)')

# Figura 1: mapa histórico (NH 2009)
plt.figure(1, figsize=(7.5, 6))
im1 = plt.pcolormesh(A2, A3, mapa_2009_NH.T, cmap=cmap_base, shading='auto')
plt.title(r'Figura 1. Topología del Espacio de Fases Histórico ($\text{NH}_{2009}$)', fontsize=11, pad=12)
plt.xlabel(r'Fase de Majorana $\alpha_2$ (rad)', fontsize=10)
plt.ylabel(r'Fase de Majorana $\alpha_3$ (rad)', fontsize=10)
cbar1 = plt.colorbar(im1, label=r'Cota Superior Asociada de $m_{eff}$ (MeV)')
plt.legend(handles=[gray_patch], loc='lower left', fontsize=9)
plt.grid(True, linestyle='--', alpha=0.3)
plt.tight_layout()
plt.savefig(FIGURAS / "03_fig1_espacio_fases_2009.png", dpi=150)

# Figura 2: jerarquías con parámetros de 2020 (NH frente a IH)
fig2, axes2 = plt.subplots(1, 2, figsize=(15, 6.2))

# Panel A: NH 2020
im2A = axes2[0].pcolormesh(A2, A3, mapa_2020_NH.T, cmap=cmap_base, shading='auto')
axes2[0].set_title(r'(a) Jerarquía Normal Moderna ($\text{NH}_{2020}$)', fontsize=11)
axes2[0].set_xlabel(r'Fase de Majorana $\alpha_2$ (rad)', fontsize=10)
axes2[0].set_ylabel(r'Fase de Majorana $\alpha_3$ (rad)', fontsize=10)
axes2[0].grid(True, linestyle='--', alpha=0.3)
fig2.colorbar(im2A, ax=axes2[0], label=r'Cota Superior de $m_{eff}$ (MeV)')
axes2[0].legend(handles=[gray_patch], loc='lower left', fontsize=9)

# Panel B: IH 2020
im2B = axes2[1].pcolormesh(A2, A3, mapa_2020_IH.T, cmap=cmap_base, shading='auto')
axes2[1].set_title(r'(b) Jerarquía Invertida Moderna ($\text{IH}_{2020}$)', fontsize=11)
axes2[1].set_xlabel(r'Fase de Majorana $\alpha_2$ (rad)', fontsize=10)
axes2[1].grid(True, linestyle='--', alpha=0.3)
fig2.colorbar(im2B, ax=axes2[1], label=r'Cota Superior de $m_{eff}$ (MeV)')
axes2[1].legend(handles=[gray_patch], loc='lower left', fontsize=9)

fig2.suptitle('Figura 2. Comparativa Topológica de las Jerarquías de Masa con Parámetros Modernos', fontsize=12, y=0.98)
plt.tight_layout()
fig2.savefig(FIGURAS / "03_fig2_jerarquias_2020.png", dpi=150)

# Figura 3: evolución 2009 a 2020 (residuos y zonas recuperadas)
fig3, axes3 = plt.subplots(1, 2, figsize=(15, 6.2))

# Panel A: Mapa de Residuos
mapa_dif_hist = mapa_2009_NH - mapa_2020_NH
max_val_hist = np.nanmax(np.abs(mapa_dif_hist)) if np.any(~np.isnan(mapa_dif_hist)) else 1.0

im3A = axes3[0].pcolormesh(A2, A3, mapa_dif_hist.T, cmap='coolwarm', vmin=-max_val_hist, vmax=max_val_hist, shading='auto')
axes3[0].set_title(r'(a) Residuo de Evolución Temporal ($\text{NH}_{2009} - \text{NH}_{2020}$)', fontsize=11)
axes3[0].set_xlabel(r'Fase de Majorana $\alpha_2$ (rad)', fontsize=10)
axes3[0].set_ylabel(r'Fase de Majorana $\alpha_3$ (rad)', fontsize=10)
axes3[0].grid(True, linestyle='--', alpha=0.3)
fig3.colorbar(im3A, ax=axes3[0], label=r'$\Delta m_{eff}$ (MeV) [Rojo: Cota 2009 > 2020 | Azul: Cota 2020 > 2009]')

# Panel B: Zonas Recuperadas
zona_iluminada = np.isnan(mapa_2009_NH) & ~np.isnan(mapa_2020_NH)
axes3[1].pcolormesh(A2, A3, zona_iluminada.T, cmap='Blues', shading='auto')
axes3[1].set_title('(b) Espacio Geométrico Recuperado (Límites antes inaccesibles)', fontsize=11)
axes3[1].set_xlabel(r'Fase de Majorana $\alpha_2$ (rad)', fontsize=10)
axes3[1].grid(True, linestyle='--', alpha=0.3)
blue_patch = mpatches.Patch(color='#3182bd', label='Zonas Iluminadas por la Precisión de 2020 (7.75%)')
axes3[1].legend(handles=[blue_patch, gray_patch], loc='lower left', fontsize=9)

fig3.suptitle('Figura 3. Impacto de la Precisión de los Parámetros de Mezcla en la Estructura del Espacio de Fases', fontsize=12, y=0.98)
plt.tight_layout()
fig3.savefig(FIGURAS / "03_fig3_evolucion_2009_2020.png", dpi=150)
