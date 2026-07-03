import math

from suspensao.modelo import SuspensaoCarro
from suspensao.graficos import plot_comparacao, plot_convergencia


def main():
    m = 350.0
    k = 20000.0
    zeta_alvo = 0.3
    c = zeta_alvo * 2.0 * math.sqrt(k * m)

    suspensao = SuspensaoCarro(m, c, k)

    print("=== Suspensao de um carro: oscilador harmonico amortecido ===")
    print(f"m = {m} kg, k = {k} N/m, c = {c:.2f} N.s/m")
    print(f"omega_n = {suspensao.omega_n:.4f} rad/s")
    print(f"zeta    = {suspensao.zeta:.4f}  ({suspensao.regime()})")
    dt_crit = suspensao.dt_critico_euler()
    print(f"dt critico de estabilidade do metodo de Euler ~= {dt_crit:.4f} s")

    x0 = 0.05
    v0 = 0.0

    plot_comparacao(
        suspensao, x0, v0, dt=0.01, tf=2.0,
        titulo="Resposta ao solavanco - passo estavel (dt = 0.01 s)",
        nome_arquivo="cenario1_dt_estavel.png",
    )

    plot_comparacao(
        suspensao, x0, v0, dt=0.15, tf=3.0,
        titulo="Resposta ao solavanco - passo instavel para Euler (dt = 0.15 s)",
        nome_arquivo="cenario2_dt_instavel.png",
    )

    lista_dt = [0.1, 0.05, 0.02, 0.01, 0.005, 0.002, 0.001]
    plot_convergencia(
        suspensao, x0, v0, t_alvo=1.0, lista_dt=lista_dt,
        nome_arquivo="convergencia_erro.png",
    )


if __name__ == "__main__":
    main()
