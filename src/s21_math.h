#ifndef SRC_S21_MATH_H_
#define SRC_S21_MATH_H_

#include <limits.h>

#define S21_EPS 1e-17
#define S21_MPI 3.141592653589793238462643383
#define S21_EPSEQ 1e-07
#define S21_LN10 2.30258509299404590109
#define S21_NAN 0.0 / 0.0
#define S21_INFINITY 1.0 / 0.0

int s21_abs(int x);
long double s21_acos(double x);
long double s21_asin(double x);
long double s21_atan(double x);
long double s21_ceil(double x);  // возвращает ближайшее целое число
long double s21_cos(double x);
long double s21_exp(double x);  // возвращает e, возведенное в заданную степень
long double s21_fabs(double x);  // вычисляет абсолютное значение числа
long double s21_floor(double x);  // возвращает ближайшее число
long double s21_fmod(double x, double y);  // остаток операции деления
long double s21_log(double x);  // вычисляет натуральный логарифм
long double s21_pow(double base,
                    double exp);  // возводит число в заданную степень
long double s21_sin(double x);
long double s21_sqrt(double x);  // вычисляет квадратный корень
long double s21_tan(double x);

//  dop func
int check_nan(double x);

#endif  // SRC_S21_MATH_H_
