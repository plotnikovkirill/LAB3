import matplotlib

try:
    matplotlib.use('TkAgg')
except:
    pass

import numpy as np
import matplotlib.pyplot as plt
from matplotlib.widgets import Slider


def simulate_nand_physics(S_val, C_val, t_switch):
    U_p = 5.0
    R = 1e3

    C_n = C_val  # Нагрузка
    C_si = C_val / 5.0  # Паразитная емкость (сделали чуть больше, чтобы видеть эффект)
    C_zi = C_val / 10.0


    t_end = 80e-9
    N = 3000
    h = t_end / N
    time = np.linspace(0, t_end, N)


    u_in = np.zeros(N)
    U_Csi1 = np.zeros(N)
    U_Csi2 = np.zeros(N)
    U_Cn = np.zeros(N)


    U_Csi1[0] = 0.0
    U_Csi2[0] = 5.0
    U_Cn[0] = 5.0

    D = 3 * C_si + 2 * C_zi + 2 * C_n
    K1 = (S_val * (C_zi + C_si + C_n)) / (C_si * D)
    K2 = (S_val * (C_zi + 2 * C_si + C_n)) / (C_si * D)
    K3 = S_val / D
    K_cross = S_val / C_si

    for n in range(N - 1):
        if time[n] >= t_switch:
            u_in[n] = 5.0
        else:
            u_in[n] = 0.0

        curr_in = u_in[n]
        uc1 = U_Csi1[n]
        uc2 = U_Csi2[n]

        if curr_in < 2.5:
            U_Cn[n + 1] = U_Cn[n] + h * ((U_p - U_Cn[n]) / (R * C_n))

            U_Csi1[n + 1] = uc1 * 0.95
            U_Csi2[n + 1] = U_Cn[n + 1] - U_Csi1[n + 1]

        else:


            if uc1 <= 0.01 and uc2 <= 0.01 and U_Cn[n] < 0.1:
                U_Csi1[n + 1] = 0;
                U_Csi2[n + 1] = 0;
                U_Cn[n + 1] = 0
                continue


            dU2_dt = (K1 * curr_in) - (K2 * curr_in) - K3 * (U_p + uc1 + uc2)
            dU1_dt = dU2_dt - (K_cross * curr_in) + (K_cross * curr_in)



            # Считаем "идеальный" делитель напряжения для текущего момента
            target_mid = U_Cn[n] / 2.0


            if uc1 < target_mid * 0.1:
                dU1_dt += (target_mid - uc1) * 1e8

            # Эйлер
            next_uc1 = uc1 + h * dU1_dt
            next_uc2 = uc2 + h * dU2_dt

            # Не уходим в минус
            U_Csi1[n + 1] = max(0, next_uc1)
            U_Csi2[n + 1] = max(0, next_uc2)

            # Выход = сумма
            U_Cn[n + 1] = U_Csi1[n + 1] + U_Csi2[n + 1]

    return time, u_in, U_Cn, U_Csi1, U_Csi2


init_S = 1.5e-3
init_C = 50e-12
init_T = 15e-9

fig, ax = plt.subplots(figsize=(10, 7))
plt.subplots_adjust(bottom=0.35)

t, u_in, u_out, u_low, u_high = simulate_nand_physics(init_S, init_C, init_T)

line_in, = ax.plot(t * 1e9, u_in, 'k--', alpha=0.3, label='Вход (Лог. 1)')
line_out, = ax.plot(t * 1e9, u_out, 'r-', lw=3, label='Выход $U_{C_н}$')
line_low, = ax.plot(t * 1e9, u_low, 'g-', lw=2, label='Нижний $U_{си1}$ (Смотри сюда!)')
line_high, = ax.plot(t * 1e9, u_high, 'b:', lw=2, label='Верхний $U_{си2}$')

ax.set_title("Лаб 2: И-НЕ с перераспределением заряда (The Bump)")
ax.set_xlabel("Время, нс")
ax.set_ylabel("Напряжение, В")
ax.grid(True)
ax.legend(loc='upper right')

ax_S = plt.axes([0.15, 0.25, 0.65, 0.03], facecolor='lightyellow')
s_S = Slider(ax_S, 'Крутизна S', 0.1, 2.0, valinit=init_S * 1000)

ax_C = plt.axes([0.15, 0.20, 0.65, 0.03], facecolor='lightyellow')
s_C = Slider(ax_C, 'Емкость C', 10, 200, valinit=init_C * 1e12)

ax_T = plt.axes([0.15, 0.15, 0.65, 0.03], facecolor='lightyellow')
s_T = Slider(ax_T, 'Старт (нс)', 5, 60, valinit=init_T * 1e9)


def update(val):
    S = s_S.val / 1000.0
    C = s_C.val * 1e-12
    T = s_T.val * 1e-9
    t_n, u_in_n, u_out_n, u_low_n, u_high_n = simulate_nand_physics(S, C, T)
    line_in.set_ydata(u_in_n)
    line_out.set_ydata(u_out_n)
    line_low.set_ydata(u_low_n)
    line_high.set_ydata(u_high_n)
    fig.canvas.draw_idle()


s_S.on_changed(update)
s_C.on_changed(update)
s_T.on_changed(update)

plt.show()