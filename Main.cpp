#include <iostream>
#include <fstream>
#include <cmath>
#include <vector>
#include <iomanip>
#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif


// Структура для зберігання вхідних даних одного варіанта
struct Variant {
    double rho, h0, beta, mu, a, S, T, Delta, Omega0, alpha, f_vyh, R, R1;
};

// Функція обчислення кутової швидкості
double Omega1(double t, double Omega0, double alpha, double T) {
    return Omega0 * (1 + alpha * cos(2 * M_PI * t / T));
}

// Функція обчислення висоти H(t)
double H_t(double t, double h0, double beta, double T) {
    if (t < T / 4.0)
        return h0 * (1 + beta);
    else if (t < 3 * T / 4.0)
        return h0;
    else
        return h0 * (1 - beta);
}

// Функція обчислення потужності N1
double N1(double t, const Variant& v, double g = 9.81) {
    double omega = Omega1(t, v.Omega0, v.alpha, v.T);
    double H = H_t(t, v.h0, v.beta, v.T);

    double base = omega * omega * v.R1 * v.R1 + 2 * g * H;

    if (base < 0 || v.a == 0 || v.S == 0) {
        std::cerr << "Error when t = " << t << ": invalid values (base < 0 or / 0)\n";
        return NAN;
    }

    double numerator = 2 * M_PI * v.rho * omega * omega * v.f_vyh * v.mu * pow(v.R, 3);
    double denominator = 3 * v.a * v.S;
    double sqrt_term = pow(base, 1.5);

    return (numerator / denominator) * sqrt_term;
}

int main() {
    std::ifstream fin("data/input.txt");  // Вхідний файл (у папці data)
    std::ofstream fout("output.txt");     // Вихідний файл

    if (!fin.is_open()) {
        std::cerr << "Error opening data/input.txt\n";
        return 1;
    }

    std::vector<Variant> variants(3);  // Масив із трьох варіантів

    // Зчитування вхідних даних
    for (int i = 0; i < 3; ++i) {
        fin >> variants[i].rho >> variants[i].h0 >> variants[i].beta >> variants[i].mu
            >> variants[i].a >> variants[i].S >> variants[i].T >> variants[i].Delta
            >> variants[i].Omega0 >> variants[i].alpha >> variants[i].f_vyh
            >> variants[i].R >> variants[i].R1;
    }

    fin.close();

    // Обчислення і виведення результатів
    for (int i = 0; i < 3; ++i) {
        fout << "=== Variant " << i + 1 << " ===\n";
        std::cout << "=== Variant " << i + 1 << " ===\n";

        for (double t = 0; t <= variants[i].T; t += variants[i].Delta) {
            double n1 = N1(t, variants[i]);

            fout << std::fixed << std::setprecision(5)
                << "t = " << t << " c, N1 = ";
            std::cout << std::fixed << std::setprecision(5)
                << "t = " << t << " c, N1 = ";

            if (std::isnan(n1)) {
                fout << "error (NaN)\n";
                std::cout << "error (NaN)\n";
            }
            else {
                fout << n1 << " watt\n";
                std::cout << n1 << " watt\n";
            }
        }

        fout << "\n";
        std::cout << "\n";
    }

    fout.close();
    return 0;
}
