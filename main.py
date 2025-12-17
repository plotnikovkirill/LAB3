import matplotlib

try:
    matplotlib.use('TkAgg')
except:
    pass

import numpy as np
import matplotlib.pyplot as plt
from matplotlib.widgets import Slider, Button



def calculate_nand_model(S_val, C_val, U_zi_val):
    """
    S_val: Крутизна (А/В)
    C_val: Емкость нагрузки (Ф) (остальные емкости масштабируем от нее)
    U_zi_val: Входное напряжение (В)
    """
    # Константы
    U_p = 5.0  # Питание

    # Емкости
    C_n = C_val  # Нагрузка
    C_zi = C_val / 10.0  # Паразитная емкость затвора
    C_si = C_val / 10.0  # Паразитная емкость стока

    # Временные параметры
    t_end = 100e-9
    N = 2000
    h = t_end / N
    time = np.linspace(0, t_end, N)

    # Массивы
    U_Csi1 = np.zeros(N)  # Нижний
    U_Csi2 = np.zeros(N)  # Верхний
    U_Cn = np.zeros(N)  # Выход

    # Начальные условия
    U_Csi1[0] = U_p / 2
    U_Csi2[0] = U_p / 2
    U_Cn[0] = U_Csi1[0] + U_Csi2[0]

    # Коэффициенты для формул
    # D - общий знаменатель системы
    D = 3 * C_si + 2 * C_zi + 2 * C_n

    # Коэффициенты K
    K1 = (S_val * (C_zi + C_si + C_n)) / (C_si * D)
    K2 = (S_val * (C_zi + 2 * C_si + C_n)) / (C_si * D)
    K3 = S_val / D
    K_cross = S_val / C_si

    # Расчет (Эйлер)
    for n in range(N - 1):
        uc1 = U_Csi1[n]
        uc2 = U_Csi2[n]

        # Если разрядились - стоп
        if uc1 <= 0 and uc2 <= 0:
            U_Csi1[n + 1] = 0
            U_Csi2[n + 1] = 0
            U_Cn[n + 1] = 0
            continue

        # Формулы из задания
        dU2_dt = (K1 * U_zi_val) - (K2 * U_zi_val) - K3 * (U_p + uc1 + uc2)
        dU1_dt = dU2_dt - (K_cross * U_zi_val) + (K_cross * U_zi_val)

        # Шаг интегрирования
        next_uc1 = uc1 + h * dU1_dt
        next_uc2 = uc2 + h * dU2_dt

        # Ограничение (не уходим в минус)
        U_Csi1[n + 1] = max(0, next_uc1)
        U_Csi2[n + 1] = max(0, next_uc2)
        U_Cn[n + 1] = U_Csi1[n + 1] + U_Csi2[n + 1]

    return time, U_Cn, U_Csi1, U_Csi2


# интерфейс

init_S = 0.3e-3  # 0.3 мА/В
init_C = 50e-12  # 50 пФ
init_Ui = 5.0  # 5 В

fig, ax = plt.subplots(figsize=(10, 7))
plt.subplots_adjust(left=0.1, bottom=0.3)  # Место под слайдеры

# Первый расчет
t, u_out, u_low, u_high = calculate_nand_model(init_S, init_C, init_Ui)

# графики
# Красный - выход
line_out, = ax.plot(t * 1e9, u_out, 'r-', lw=2, label=r'$U_{вых} (Нагрузка)$')
# Зеленый - нижний транзистор (ТОЛСТЫЙ)
line_low, = ax.plot(t * 1e9, u_low, 'g-', lw=5, alpha=0.4, label=r'$U_{си1} (Нижний)$')
# Синий - верхний транзистор (Пунктир внутри зеленого)
line_high, = ax.plot(t * 1e9, u_high, 'b--', lw=1.5, label=r'$U_{си2} (Верхний)$')

ax.set_title("Лаб 2: Интерактивная модель И-НЕ")
ax.set_xlabel("Время, нс")
ax.set_ylabel("Напряжение, В")
ax.set_ylim(-0.5, 6)
ax.set_xlim(0, 100)
ax.grid(True)
ax.legend()


# Слайдеры

axcolor = 'lightgoldenrodyellow'

# Слайдер S (Крутизна)
ax_S = plt.axes([0.15, 0.20, 0.65, 0.03], facecolor=axcolor)
s_S = Slider(ax_S, 'Крутизна S [мА/В]', 0.05, 5.0, valinit=init_S * 1000)

# Слайдер C (Емкость)
ax_C = plt.axes([0.15, 0.15, 0.65, 0.03], facecolor=axcolor)
s_C = Slider(ax_C, 'Емкость Cн [пФ]', 10, 200, valinit=init_C * 1e12)

# Слайдер U (Вход)
ax_U = plt.axes([0.15, 0.10, 0.65, 0.03], facecolor=axcolor)
s_U = Slider(ax_U, 'Вход U_зи [В]', 0.0, 5.0, valinit=init_Ui)



# update

def update(val):
    S_new = s_S.val / 1000.0
    C_new = s_C.val * 1e-12
    U_new = s_U.val

    # Пересчитываем
    t_new, u_out_new, u_low_new, u_high_new = calculate_nand_model(S_new, C_new, U_new)

    # Обновляем линии
    line_out.set_ydata(u_out_new)
    line_low.set_ydata(u_low_new)
    line_high.set_ydata(u_high_new)

    fig.canvas.draw_idle()


s_S.on_changed(update)
s_C.on_changed(update)
s_U.on_changed(update)

plt.show()