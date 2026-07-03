import math


class SuspensaoCarro:
    def __init__(self, m, c, k):
        self.m = m
        self.c = c
        self.k = k
        self.omega_n = math.sqrt(k / m)
        self.zeta = c / (2.0 * math.sqrt(k * m))

    def derivadas(self, x, v):
        return v, -(self.c * v + self.k * x) / self.m

    def regime(self):
        if self.zeta < 1.0:
            return "subamortecido"
        if abs(self.zeta - 1.0) < 1e-9:
            return "criticamente amortecido"
        return "superamortecido"

    def solucao_analitica(self, t, x0, v0):
        wn, zeta = self.omega_n, self.zeta
        if zeta < 1.0:
            wd = wn * math.sqrt(1.0 - zeta ** 2)
            A = x0
            B = (v0 + zeta * wn * x0) / wd
            return math.exp(-zeta * wn * t) * (A * math.cos(wd * t) + B * math.sin(wd * t))
        if abs(zeta - 1.0) < 1e-9:
            C1 = x0
            C2 = v0 + wn * x0
            return (C1 + C2 * t) * math.exp(-wn * t)
        wd = wn * math.sqrt(zeta ** 2 - 1.0)
        r1 = -wn * zeta + wd
        r2 = -wn * zeta - wd
        C1 = (v0 - r2 * x0) / (r1 - r2)
        C2 = x0 - C1
        return C1 * math.exp(r1 * t) + C2 * math.exp(r2 * t)

    def dt_critico_euler(self):
        a = self.zeta * self.omega_n
        b = self.omega_n * math.sqrt(abs(1.0 - self.zeta ** 2))
        return 2.0 * a / (a ** 2 + b ** 2)
