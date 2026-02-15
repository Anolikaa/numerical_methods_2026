import numpy as np


def f(x):
    return np.sin(x)


def progonka(y, h, N):
    # В Python індексація починається з 0, тому адаптуємо розміри
    alfa = np.zeros(N + 1)
    beta = np.zeros(N + 1)
    hamma = np.zeros(N + 1)
    delta = np.zeros(N + 1)
    A = np.zeros(N + 1)
    B = np.zeros(N + 1)
    c = np.zeros(N + 1)

    beta[1] = 1.0

    # Прямий хід: обчислення коефіцієнтів
    for i in range(2, N + 1):
        alfa[i] = h[i - 1]
        beta[i] = 2 * (h[i - 1] + h[i])
        hamma[i] = h[i]
        delta[i] = 3 * (((y[i] - y[i - 1]) / h[i]) - ((y[i - 1] - y[i - 2]) / h[i - 1]))

    hamma[N] = 0.0
    A[1] = -hamma[1] / beta[1]
    B[1] = delta[1] / beta[1]

    for i in range(2, N):
        denom = alfa[i] * A[i - 1] + beta[i]
        A[i] = -hamma[i] / denom
        B[i] = (delta[i] - alfa[i] * B[i - 1]) / denom

    # Зворотний хід: знаходження невідомих c[i]
    c[N] = (delta[N] - alfa[N] * B[N - 1]) / (alfa[N] * A[N - 1] + beta[N])
    for i in range(N, 1, -1):
        c[i - 1] = A[i - 1] * c[i] + B[i - 1]

    return c


def main():
    # Зчитування даних за допомогою numpy
    try:
        # Припускаємо формат: індекс, x, y, h
        data = np.loadtxt("input.txt")
        x = data[:, 1]
        y = data[:, 2]
        h = data[:, 3]
        N = len(x) - 1
    except Exception as e:
        print(f"Помилка зчитування файлу: {e}")
        return

    x0 = x[0]
    hh = h[0]
    c = progonka(y, h, N)

    # Розрахунок коефіцієнтів сплайна
    a = np.zeros(N + 1)
    b = np.zeros(N + 1)
    d = np.zeros(N + 1)

    for i in range(1, N):
        a[i] = y[i - 1]
        b[i] = (y[i] - y[i - 1]) / h[i] - (h[i] / 3) * (c[i + 1] + 2 * c[i])
        d[i] = (c[i + 1] - c[i]) / (3 * h[i])

    a[N] = y[N - 1]
    b[N] = (y[N] - y[N - 1]) / h[N] - (2.0 / 3.0) * h[N] * c[N]
    d[N] = -c[N] / (3 * h[N])

    # Генерація точок для виводу (20 точок на інтервал)
    with open("output.txt", "w") as f_out:
        step = hh / 20.0
        j = 1
        for i in range(20 * (N - 1) + 1):
            xm = x0 + i * step
            ym = f(xm)

            dx = xm - x[j - 1]
            # Формула кубічного сплайна: S(x) = a + b(x-x0) + c(x-x0)^2 + d(x-x0)^3
            s = a[j] + b[j] * dx + c[j] * (dx ** 2) + d[j] * (dx ** 3)
            eps = abs(s - ym)

            f_out.write(f"[{i},{j}]\t{xm:.6e}\t{ym:.6e}\t{s:.6e}\t{eps:.6e}\n")

            if i != 0 and i % 20 == 0:
                j += 1

    print("Готово! Результати збережено в output.txt")


if __name__ == "__main__":
    main()
