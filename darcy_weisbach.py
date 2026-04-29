import math
import numpy as np
from friccion import factor_friccion

# ---- DATOS ----
C = 3000 # l/hr
Longitud_recta = 21.68 # m
Altura = 10.28 # m
velocidad_deseada = 0.7 # m/s

# ---- Accesorios ----
Codos90 = 7
Te_recta= 0
Te_lateral = 0
LLP = 1
Valvula_retencion = 1

# ---- CALCULO ----

Q = C / (1000*3600) # m^3/s
eps = 1.5e-6 
f=0.025

def v_med(Q, D):
    v_med = (4*Q) / (math.pi * D**2)
    return v_med

def reynolds(Q, D):
    Re=(v_med(Q, D))*D*(1/1.0038e-6)
    return Re

def perdida_carga(f, D, v):
    R = (f * (v**2)) / (2*D*9.81)
    return R

def l_equivalente(ki,D,f):
    L_eq = (ki * D) / f
    return L_eq

def Hf(h, L, Leq, j):
    Hf = h + (L + Leq) * j
    return Hf

# Te recta, Te lateral, Codo 90°, LLP, V-Retencion
K = [1.5, 2.5, 1.5, 0.5, 5] 
D = [float(x) for x in open("Calculo/lista-diametros.txt", "r").readlines()]
Re = [reynolds(Q, D[i]) for i in range(len(D))]
v = [v_med(Q, D[i]) for i in range(len(D))]
deltav = [abs(v[i] - velocidad_deseada) for i in range(len(D))]

print("---------------------------------------------------------------------")
print("     ","Longitudes equivalentes de accesorios en metros")
print("Diametro", " ", "Te recta", " ", "Te lateral", "  ", "Codo 90°", "   ", "LLP" ," ", "V-Retencion")
for j in range(len(D)):
    fr = [factor_friccion(f, D[j], eps, Re[j]) for j in range(len(D))]
    L_eq0 = l_equivalente(K[0], D[j], fr[j])
    L_eq1 = l_equivalente(K[1], D[j], fr[j])
    L_eq2 = l_equivalente(K[2], D[j], fr[j])
    L_eq3 = l_equivalente(K[3], D[j], fr[j])
    L_eq4 = l_equivalente(K[4], D[j], fr[j])
    print(f"{D[j]:.3f}", "     ", f"{L_eq0:.3f} m","   ", f"{L_eq1:.3f} m", "   ", f"{L_eq2:.3f} m", "  ", f"{L_eq3:.3f} m", "  ", f"{L_eq4:.3f} m")

print("---------------------------------------------------------------------")
print("diam", " ", "Velocidad [m/s]", " ", "Hf [m]", "", "Perdida [m/m]", "", "L-eq [m]")
for i in range(len(D)):
    if deltav[i] == np.min(deltav) and v[i] > 0.5:
        fr = [factor_friccion(f, D[i], eps, Re[i]) for i in range(len(D))]
        R = [perdida_carga(fr[i], D[i], v[i]) for i in range(len(D))]
        L_eq0 = l_equivalente(K[0], D[i], fr[i])
        L_eq1 = l_equivalente(K[1], D[i], fr[i])
        L_eq2 = l_equivalente(K[2], D[i], fr[i])
        L_eq3 = l_equivalente(K[3], D[i], fr[i])
        L_eq4 = l_equivalente(K[4], D[i], fr[i])
        L_eqtotal = (L_eq0*Te_recta) + (L_eq1*Te_lateral) + (L_eq2*Codos90) + (L_eq3*LLP) + (L_eq4*Valvula_retencion)
        H = Hf(Altura, Longitud_recta, L_eqtotal, R[i])
        HP = (H*Q*1000*0.013)/0.4
        print(f"{D[i]}","    ", f"{v[i]:.3f}","       ", f"{H:.3f}", "    ", "R =", f"{R[i]:.3f}", "  ", f"{L_eqtotal:.3f}")
        print("Bomba »", f"{HP:.3f} HP")
    elif deltav[i] == np.min(deltav):
        print("No se encontró un diámetro con velocidad mayor a 0.5 m/s.")
        fr = [factor_friccion(f, D[i], eps, Re[i]) for i in range(len(D))]
        R = [perdida_carga(fr[i], D[i], v[i]) for i in range(len(D))]
        L_eq0 = l_equivalente(K[0], D[i-1], fr[i-1])
        L_eq1 = l_equivalente(K[1], D[i-1], fr[i-1])
        L_eq2 = l_equivalente(K[2], D[i-1], fr[i-1])
        L_eq3 = l_equivalente(K[3], D[i-1], fr[i-1])
        L_eq4 = l_equivalente(K[4], D[i], fr[i])
        L_eqtotal = (L_eq0*Te_recta) + (L_eq1*Te_lateral) + (L_eq2*Codos90) + (L_eq3*LLP) + (L_eq4*Valvula_retencion)
        H = Hf(Altura, Longitud_recta, L_eqtotal, R[i-1])
        HP = (H*Q*1000*0.013)/0.4
        print(f"{D[i-1]}","   ", f"{v[i-1]:.3f}","  ", f"{H:.3f}", "    ", "R =", f"{R[i-1]:.3f}", "    ", f"{L_eqtotal:.3f}")
        print("Bomba »", f"{HP:.3f} HP")
print("---------------------------------------------------------------------")
