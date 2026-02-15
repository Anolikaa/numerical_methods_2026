import numpy as np


def f(x):
    return np.sin(x)


def generate_input():
    N = 50
    x0 = 0.0
    xn = 1.0

    # Створюємо масив x з N+1 точкою від x0 до xn
    x = np.linspace(x0, xn, N + 1)
    hh = (xn - x0) / N

    # Відкриваємо файл для запису
    with open("input.txt", "w") as finput:
        for i in range(N + 1):
            xi = x[i]
            yi = f(xi)
            # Формат: індекс, x, y, h
            finput.write(f"{i}\t{xi:.6e}\t{yi:.6e}\t{hh:.6e}\n")

    print("Файл input.txt успішно згенеровано! (DONE)")


if __name__ == "__main__":
    generate_input()