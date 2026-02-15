import numpy as np
import matplotlib.pyplot as plt


def f(x):
    return np.sin(x)


def progonka(y, h, N):
    alfa = np.zeros(N + 1)
    beta = np.zeros(N + 1)
    hamma = np.zeros(N + 1)
    delta = np.zeros(N + 1)
    A = np.zeros(N + 1)
    B = np.zeros(N + 1)
    c = np.zeros(N + 1)

    beta[1] = 1.0
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

    c[N] = (delta[N] - alfa[N] * B[N - 1]) / (alfa[N] * A[N - 1] + beta[N])
    for i in range(N, 1, -1):
        c[i - 1] = A[i - 1] * c[i] + B[i - 1]

    return c


def main():
    try:
        data = np.loadtxt("input.txt")
        x_nodes = data[:, 1]
        y_nodes = data[:, 2]
        h = data[:, 3]
        N = len(x_nodes) - 1
    except Exception as e:
        print(f"Помилка: {e}. Перевірте файл input.txt")
        return

    c = progonka(y_nodes, h, N)

    a = np.zeros(N + 1)
    b = np.zeros(N + 1)
    d = np.zeros(N + 1)

    for i in range(1, N):
        a[i] = y_nodes[i - 1]
        b[i] = (y_nodes[i] - y_nodes[i - 1]) / h[i] - (h[i] / 3) * (c[i + 1] + 2 * c[i])
        d[i] = (c[i + 1] - c[i]) / (3 * h[i])

    a[N] = y_nodes[N - 1]
    b[N] = (y_nodes[N] - y_nodes[N - 1]) / h[N] - (2.0 / 3.0) * h[N] * c[N]
    d[N] = -c[N] / (3 * h[N])

    # Списки для графіків
    x_spline = []
    y_spline = []

    step = h[0] / 20.0
    j = 1
    for i in range(20 * (N - 1) + 1):
        xm = x_nodes[0] + i * step
        dx = xm - x_nodes[j - 1]
        s = a[j] + b[j] * dx + c[j] * (dx ** 2) + d[j] * (dx ** 3)

        x_spline.append(xm)
        y_spline.append(s)

        if i != 0 and i % 20 == 0:
            j += 1

    # ВІЗУАЛІЗАЦІЯ
    plt.figure(figsize=(10, 6))
    plt.plot(x_spline, y_spline, label='Кубічний сплайн', color='blue', linewidth=2)
    plt.scatter(x_nodes, y_nodes, color='red', label='Вузли (дані)', zorder=5)

    # Справжній синус для порівняння
    x_real = np.linspace(min(x_nodes), max(x_nodes), 100)
    plt.plot(x_real, np.sin(x_real), '--', label='sin(x)', color='green', alpha=0.5)

    plt.title('Інтерполяція кубічним сплайном')
    plt.xlabel('x')
    plt.ylabel('y')
    plt.legend()
    plt.grid(True)

    # Збереження графіка та показ
    plt.savefig("spline_plot.png")
    print("Графік збережено як 'spline_plot.png'")
    plt.show()


if __name__ == "__main__":
    main()