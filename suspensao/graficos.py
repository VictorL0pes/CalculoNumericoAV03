import math
import os

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

from suspensao.metodos_numericos import metodo_euler, metodo_rk4, erro_absoluto

PASTA_SAIDA = os.path.join(
    os.path.dirname(os.path.dirname(os.path.abspath(__file__))), "graficos"
)


def salvar(fig, nome):
    os.makedirs(PASTA_SAIDA, exist_ok=True)
    caminho = os.path.join(PASTA_SAIDA, nome)
    fig.savefig(caminho, dpi=150, bbox_inches="tight")
    print(f"Grafico salvo em: {caminho}")


def plot_comparacao(sistema, x0, v0, dt, tf, titulo, nome_arquivo):
    ts_e, xs_e, _ = metodo_euler(sistema, x0, v0, 0.0, tf, dt)
    ts_r, xs_r, _ = metodo_rk4(sistema, x0, v0, 0.0, tf, dt)
    xs_a = [sistema.solucao_analitica(t, x0, v0) for t in ts_e]

    fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(9, 8), sharex=True)

    ax1.plot(ts_e, xs_a, "k-", linewidth=2, label="Solucao analitica")
    ax1.plot(ts_e, xs_e, "r--", linewidth=1.5, label=f"Euler (dt={dt}s)")
    ax1.plot(ts_r, xs_r, "b-.", linewidth=1.5, label=f"RK4 (dt={dt}s)")
    ax1.set_ylabel("Deslocamento x(t) [m]")
    ax1.set_title(titulo)
    ax1.legend()
    ax1.grid(True, alpha=0.3)

    erro_e = erro_absoluto(ts_e, xs_e, sistema, x0, v0)
    erro_r = erro_absoluto(ts_r, xs_r, sistema, x0, v0)
    ax2.semilogy(ts_e, [max(e, 1e-16) for e in erro_e], "r--", label="Erro absoluto - Euler")
    ax2.semilogy(ts_r, [max(e, 1e-16) for e in erro_r], "b-.", label="Erro absoluto - RK4")
    ax2.set_xlabel("Tempo [s]")
    ax2.set_ylabel("Erro absoluto |x_num - x_exato| [m]")
    ax2.legend()
    ax2.grid(True, alpha=0.3)

    salvar(fig, nome_arquivo)
    plt.close(fig)


def plot_convergencia(sistema, x0, v0, t_alvo, lista_dt, nome_arquivo):
    erros_euler, erros_rk4 = [], []
    for dt in lista_dt:
        _, xs_e, _ = metodo_euler(sistema, x0, v0, 0.0, t_alvo, dt)
        _, xs_r, _ = metodo_rk4(sistema, x0, v0, 0.0, t_alvo, dt)
        x_exato = sistema.solucao_analitica(t_alvo, x0, v0)
        erros_euler.append(abs(xs_e[-1] - x_exato))
        erros_rk4.append(abs(xs_r[-1] - x_exato))

    print("\ndt        erro Euler      erro RK4        ordem Euler   ordem RK4")
    for i, dt in enumerate(lista_dt):
        if i == 0:
            print(f"{dt:<9} {erros_euler[i]:<15.6e} {erros_rk4[i]:<15.6e} {'-':<13} -")
        else:
            razao = math.log(lista_dt[i - 1] / dt)
            ord_e = math.log(erros_euler[i - 1] / max(erros_euler[i], 1e-16)) / razao
            ord_r = math.log(erros_rk4[i - 1] / max(erros_rk4[i], 1e-16)) / razao
            print(f"{dt:<9} {erros_euler[i]:<15.6e} {erros_rk4[i]:<15.6e} {ord_e:<13.2f} {ord_r:.2f}")

    fig, ax = plt.subplots(figsize=(7, 6))
    ax.loglog(lista_dt, erros_euler, "ro-", label="Euler (ordem teorica 1)")
    ax.loglog(lista_dt, erros_rk4, "bs-", label="RK4 (ordem teorica 4)")
    ax.set_xlabel("Passo de integracao dt [s]")
    ax.set_ylabel(f"Erro absoluto em t={t_alvo}s [m]")
    ax.set_title("Convergencia dos metodos numericos")
    ax.legend()
    ax.grid(True, which="both", alpha=0.3)
    salvar(fig, nome_arquivo)
    plt.close(fig)
