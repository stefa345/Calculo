import math
tol = 2e-3            # tolerancia 
tol_g = 2e-3          # tolerancia 
max_iter = 100

def g(f, D, eps, Re):
    """g(f) de Colebrook. Devuelve +inf si f<=0 para evitar errores."""
    if f <= 0:
        return float('inf')
    return 1.0/math.sqrt(f) + 2.0*math.log10((eps/D)/3.7 + 2.51/(Re*math.sqrt(f)))

def g_prime(f, D, eps, Re):
    """Derivada analítica g'(f)."""
    if f <= 0:
        return float('inf')
    # S(f) = eps/D/3.7 + 2.51/(Re*sqrt(f))
    S = (eps / D) / 3.7 + 2.51 / (Re * math.sqrt(f))
    # S'(f) = derivative of 2.51/(Re*sqrt(f)) = - (2.51)/(2*Re) * f^(-3/2)
    S_prime = - (2.51) / (2.0 * Re) * f**(-1.5)
    # g'(f) = -1/2 f^{-3/2} + 2/(ln(10)) * S'(f)/S
    term1 = -0.5 * f**(-1.5)
    term2 = (2.0 / math.log(10.0)) * (S_prime / S)
    return term1 + term2

def factor_friccion(f, D, eps, Re):
    for k in range(1, max_iter+1):
        gf = g(f, D, eps, Re)
        gpf = g_prime(f, D, eps, Re)
        if gpf == 0 or math.isinf(gpf):
            raise ZeroDivisionError(f"Derivada nula o inválida en iter {k}: g'(f)={gpf}")

        f_next = f - gf / gpf

        # asegurar f_next positivo (si da negativo, ajustar paso)
        if f_next <= 0:
            # reducir paso (backtracking simple)
            factor = 0.5
            while f_next <= 0 and factor > 1e-6:
                f_next = f - factor * (gf / gpf)
                factor *= 0.5
            if f_next <= 0:
                raise ValueError(f"No se pudo encontrar un paso positivo en iter {k}")

        # calcular error relativo y verificar tolerancia
        err = abs(f_next - f) / abs(f) if f != 0 else None

        if err is not None and err < tol:
            break

        f = f_next
    else:
        print("\nNo se alcanzó convergencia en el número máximo de iteraciones.")
    return f