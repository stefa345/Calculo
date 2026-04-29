import math
import numpy as np
from friccion import factor_friccion

# -------- DATOS --------
C = 3000            # l/hr
Longitud_recta = 21.68   # m
Altura = 10.28      # m
velocidad_deseada = 3  # m/s
velocidad_min = 0.5      # m/s

eps = 1.5e-6        # rugosidad absoluta en metros
f0 = 0.025          # fricción inicial

# -------- Accesorios --------
accesorios = {
    "Te_recta": 0,
    "Te_lateral": 0,
    "Codo_90": 7,
    "LLP": 1,
    "Valvula_retencion": 1
}

K = {
    "Te_recta": 1.5,
    "Te_lateral": 2.5,
    "Codo_90": 1.5,
    "LLP": 0.5,
    "Valvula_retencion": 5
}

# -------- FUNCIONES --------
Q = C / (1000*3600) # m^3/s
eps = 1.5e-6 
f=0.025

def caudal(C):
    return C / (1000*3600)

def v_med(Q, D):
    v_med = (4*Q) / (math.pi * D**2)
    return v_med

def reynolds(Q, D, nu=1.0038e-6):
    Re=(v_med(Q, D))*D*(1/nu)
    return Re

def perdida_carga(f, D, v):
    R = (f * (v**2)) / (2*D*9.81)
    return R

def l_equivalente(D,f):
    total = 0
    for key in K:
        Leq = (K[key] * D) / f
        total += Leq * accesorios[key]
    return total

def Hf(h, L, Leq, j):
    Hf = h + (L + Leq) * j
    return Hf

def pot_bomba(H,Q, eficiencia=0.4):
    HP = (H*Q*1000*0.013)/eficiencia
    return HP

# -------- CALCULO --------
Q = caudal(C)

with open("lista-diametros.txt", "r") as f:
    D = np.array([float(x) for x in f.readlines()])

v = np.array([v_med(Q,d) for d in D])
Re = np.array([reynolds(v[i],D[i]) for i in range(len(D))])
fr = np.array([factor_friccion(f0, D[i], eps, Re[i]) for i in range(len(D))])
j = np.array([perdida_carga(fr[i],D[i],v[i]) for i in range(len(D))])

deltav = np.abs(v - velocidad_deseada)

verificar = np.where(v >= velocidad_min)[0]
if len(verificar) == 0:
    print("No existe diámetro con velocidad mayor a la mínima")
    idx = np.argmin(deltav)
else:
    idx = verificar[np.argmin(deltav[verificar])]

# -------- RESULTADOS --------
D_sel = D[idx]
v_sel = v[idx]
f_sel = fr[idx]
j_sel = j[idx]

Leq_total = l_equivalente(D_sel, f_sel)
H = Hf(Altura, Longitud_recta, Leq_total, j_sel)
HP = pot_bomba(H, Q)

# -------- SALIDA --------
print("---------------------------------------------------------------------")
print("          ","Longitudes equivalentes de accesorios en metros")
print("Diametro", " ", "Te recta", "  ", "Te lateral", "   ", "Codo 90°", "    ", "LLP" ,"   ", "V-Retencion")
for x in range(len(D)):
    leq0 = (K["Te_recta"] * D[x]) / fr[x]
    leq1 = (K["Te_lateral"] * D[x]) / fr[x]
    leq2 = (K["Codo_90"] * D[x]) / fr[x]
    leq3 = (K["LLP"] * D[x]) / fr[x]
    leq4 = (K["Valvula_retencion"] * D[x]) / fr[x]
    print(f"{D[x]:.3f}", "     ", f"{leq0:.3f}", "      ", f"{leq1:.3f}", "        ", f"{leq2:.3f}", "      ", f"{leq3:.3f}", "    ", f"{leq4:.3f}")


print("-----------------------------------")
print(f"Diámetro: {D_sel:.3f} m")
print(f"Velocidad: {v_sel:.3f} m/s")
print(f"Pérdida: {j_sel:.3f} m/m")
print(f"Longitud equivalente: {Leq_total:.3f} m")
print(f"Hf: {H:.3f} m")
print(f"Potencia de bomba: {HP:.3f} HP")
print("-----------------------------------")
