# Reproduccion de la cota de la masa efectiva del neutrino muonico de Majorana
# a partir del formalismo de Delepine et al. (2010), con los parametros de 2009.
# Autor: Cesar Elias Pena Aparicio
#
# Uso:     python scripts/01_reproduccion_cota_2010.py
# Salida:  figuras/01_espacio_fases_cota_2010.png y resumen por consola.

from pathlib import Path

import numpy as np
import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches

FIGURAS = Path(__file__).resolve().parents[1] / "figuras"
FIGURAS.mkdir(exist_ok=True)

# Paso 1: parametros de 2009, los usados en Delepine et al. (2010)

# Energia del haz de MINOS
E_nu_MINOS_GeV = 2.0             # GeV
E_nu = E_nu_MINOS_GeV * 1e9      # eV

# Diferencias de masas cuadradas
dm21_sq_2009 = 7.6e-5            # eV^2
dm32_sq_2009 = 2.43e-3           # eV^2

# Angulos de mezcla
theta12_2009 = np.arcsin(np.sqrt(0.304))   # sin^2 theta12 = 0.304
theta23_2009 = np.pi / 4                   # mezcla maxima, sin^2 theta23 = 0.5
theta13_2009 = 0.0                         # theta13 = 0 en el analisis de Delepine

# Fase de Dirac: irrelevante porque theta13 = 0
delta_cp_2009 = 0.0

# Linea de base de MINOS
L_MINOS = 735.0                  # km

print("Reproduccion de la cota de 2010")

# Paso 2: matriz PMNS y fases de oscilacion precalculadas
c12, s12 = np.cos(theta12_2009), np.sin(theta12_2009)
c13, s13 = np.cos(theta13_2009), np.sin(theta13_2009)
c23, s23 = np.cos(theta23_2009), np.sin(theta23_2009)

U_std = np.array([
    [c12*c13, s12*c13, s13*np.exp(-1j*delta_cp_2009)],
    [-s12*c23 - c12*s23*s13*np.exp(1j*delta_cp_2009), c12*c23 - s12*s23*s13*np.exp(1j*delta_cp_2009), s23*c13],
    [ s12*s23 - c12*c23*s13*np.exp(1j*delta_cp_2009), -c12*s23 - s12*c23*s13*np.exp(1j*delta_cp_2009), c23*c13]
])

DM2 = np.zeros((3, 3))
DM2[0, 1] = -dm21_sq_2009
DM2[0, 2] = -(dm21_sq_2009 + dm32_sq_2009)
DM2[1, 2] = -dm32_sq_2009

# Fase de oscilacion: Delta m^2 L / 4E = 1.267 * Dm2[eV^2] * L[km] / E[GeV]
fase_base = np.zeros((3, 3))
for i in range(3):
    for j in range(i+1, 3):
        fase_base[i, j] = 1.267 * DM2[i, j] * L_MINOS / E_nu_MINOS_GeV

# Paso 3: barrido de las fases de Majorana y de la masa mas ligera
res = 50                                  # resolucion de la malla de fases
alpha_vals = np.linspace(0, 2*np.pi, res)
m_light_vals = np.logspace(6, 8, 250)     # masas de 1 MeV a 100 MeV, en eV

# Cota de m_eff (MeV) para cada par de fases
matriz_cotas = np.zeros((res, res))

print("Barrido del espacio de fases...")

for i_a2, a2 in enumerate(alpha_vals):
    for i_a3, a3 in enumerate(alpha_vals):
        P = np.diag([1.0, np.exp(1j*a2), np.exp(1j*a3)])
        U = U_std @ P

        cota_alcanzada = np.nan   # punto ciego por omision

        # Se aumenta la masa hasta que se viola el limite experimental
        for ml in m_light_vals:
            masses = np.array([ml, np.sqrt(ml**2 + dm21_sq_2009), np.sqrt(ml**2 + dm21_sq_2009 + dm32_sq_2009)])

            U_mu = U[1, :]
            m_eff_eV = np.abs(np.sum(U_mu**2 * masses))

            osc_term = 0.0
            for i in range(3):
                for j in range(i+1, 3):
                    U2i_U2jstar = U[1, i]**2 * np.conj(U[1, j])**2
                    sin2_arg = np.sin(fase_base[i, j])**2
                    sin_arg = np.sin(2.0 * fase_base[i, j])

                    osc_term += 4 * np.real(U2i_U2jstar) * masses[i] * masses[j] * sin2_arg
                    osc_term += 2 * np.imag(U2i_U2jstar) * masses[i] * masses[j] * sin_arg

            rhs = m_eff_eV**2 - osc_term
            lhs = 0.001 * E_nu**2

            if rhs > lhs:
                # Limite violado: la masa de este paso fija la cota
                cota_alcanzada = m_eff_eV / 1e6   # MeV
                break

        matriz_cotas[i_a2, i_a3] = cota_alcanzada

# Paso 4: estadisticas sobre los puntos que no son ciegos
cotas_validas = matriz_cotas[~np.isnan(matriz_cotas)]

if len(cotas_validas) > 0:
    print("\nResultados del espacio de fases")
    print(f"Cota mas estricta (minima):        {np.min(cotas_validas):.2f} MeV")
    print(f"Cota mediana:                      {np.median(cotas_validas):.2f} MeV")
    print(f"Cota conservadora (maxima):        {np.max(cotas_validas):.2f} MeV")
else:
    print("Todo el espacio de fases es divergente en este rango.")

# Paso 5: mapa de la cota sobre el plano de fases de Majorana
plt.figure(figsize=(8, 6))
cmap = plt.cm.viridis.copy()
cmap.set_bad(color="gray")

A2, A3 = np.meshgrid(alpha_vals, alpha_vals)
im = plt.pcolormesh(A2, A3, matriz_cotas.T, cmap=cmap, shading="auto")

plt.colorbar(im, label="Cota maxima de $m_{eff}$ (MeV)")
plt.title("Cota de la masa efectiva segun las fases de Majorana (parametros de 2009)")
plt.xlabel("Fase $\\alpha_2$ (rad)")
plt.ylabel("Fase $\\alpha_3$ (rad)")

gray_patch = mpatches.Patch(color="gray", label="Puntos ciegos (divergencia)")
plt.legend(handles=[gray_patch], loc="upper right")

plt.tight_layout()
plt.savefig(FIGURAS / "01_espacio_fases_cota_2010.png", dpi=150)
plt.close()

# Paso 6: metodo directo, independiente de las masas
# Con masas pequenas, la cota conservadora es la raiz de 0.001 * E_nu^2
cota_conservadora_MeV = np.sqrt(0.001) * E_nu / 1e6
print("\nMetodo directo, independiente de las masas")
print(f"Cota conservadora estimada:        {cota_conservadora_MeV:.1f} MeV")

# Paso 7: aproximacion analitica de Delepine (m1 ~ m2 y coeficiente A > 0)
print("\nAproximacion analitica de Delepine (m1 ~ m2, A > 0)")

cotas_aprox = []

for i_a2, a2 in enumerate(alpha_vals):
    for i_a3, a3 in enumerate(alpha_vals):
        P = np.diag([1.0, np.exp(1j*a2), np.exp(1j*a3)])
        U = U_std @ P
        U_mu = U[1, :]

        # Masas degeneradas: m^2 se factoriza
        termino_masa = np.abs(np.sum(U_mu**2))**2

        # Terminos oscilatorios que dependen solo de las fases
        termino_osc = 0.0
        for i in range(3):
            for j in range(i+1, 3):
                U2i_U2jstar = U_mu[i]**2 * np.conj(U_mu[j])**2
                sin2_arg = np.sin(fase_base[i, j])**2
                sin_arg = np.sin(2.0 * fase_base[i, j])

                termino_osc += 4 * np.real(U2i_U2jstar) * sin2_arg
                termino_osc += 2 * np.imag(U2i_U2jstar) * sin_arg

        # Coeficiente geometrico A que multiplica a m^2
        coeficiente_total = termino_masa - termino_osc

        # Solo A > 0: se descarta la interferencia destructiva total
        if coeficiente_total > 0.0:
            cota_m_eV = np.sqrt((0.001 * E_nu**2) / coeficiente_total)
            cotas_aprox.append(cota_m_eV / 1e6)

cotas_aprox = np.array(cotas_aprox)

# Cota analitica sobre la masa efectiva <m_mumu>, sin cancelaciones
cota_m_efectiva_MeV = np.sqrt(0.001) * E_nu / 1e6

print(f"Cota analitica sobre <m_mumu>:     {cota_m_efectiva_MeV:.2f} MeV  (Delepine et al. 2010: ~64 MeV)")
print(f"Cota minima sobre la masa:         {np.min(cotas_aprox):.2f} MeV")
print(f"Cota mediana sobre la masa:        {np.median(cotas_aprox):.2f} MeV")
print(f"Cota maxima sobre la masa:         {np.max(cotas_aprox):.2f} MeV")
