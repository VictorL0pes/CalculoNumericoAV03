def metodo_euler(sistema, x0, v0, t0, tf, dt):
    n = int(round((tf - t0) / dt))
    ts, xs, vs = [t0], [x0], [v0]
    x, v, t = x0, v0, t0
    for _ in range(n):
        dx, dv = sistema.derivadas(x, v)
        x = x + dt * dx
        v = v + dt * dv
        t = t + dt
        ts.append(t)
        xs.append(x)
        vs.append(v)
    return ts, xs, vs


def metodo_rk4(sistema, x0, v0, t0, tf, dt):
    n = int(round((tf - t0) / dt))
    ts, xs, vs = [t0], [x0], [v0]
    x, v, t = x0, v0, t0
    for _ in range(n):
        k1x, k1v = sistema.derivadas(x, v)
        k2x, k2v = sistema.derivadas(x + 0.5 * dt * k1x, v + 0.5 * dt * k1v)
        k3x, k3v = sistema.derivadas(x + 0.5 * dt * k2x, v + 0.5 * dt * k2v)
        k4x, k4v = sistema.derivadas(x + dt * k3x, v + dt * k3v)
        x = x + (dt / 6.0) * (k1x + 2 * k2x + 2 * k3x + k4x)
        v = v + (dt / 6.0) * (k1v + 2 * k2v + 2 * k3v + k4v)
        t = t + dt
        ts.append(t)
        xs.append(x)
        vs.append(v)
    return ts, xs, vs


def erro_absoluto(ts, xs, sistema, x0, v0):
    return [abs(x - sistema.solucao_analitica(t, x0, v0)) for t, x in zip(ts, xs)]
