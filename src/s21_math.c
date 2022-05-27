#include "s21_math.h"

int s21_abs(int x) {  // 1
  if (x < 0) {
    x *= -1;
  }
  return x;
}

long double s21_acos(double x) {  // 2
  long double res;
  if (x < -1 || x > 1 || check_nan(x))
    res = S21_NAN;
  else if (s21_fabs(x - 1.) < S21_EPS)
    res = 0;
  else if (s21_fabs(x + 1.) < S21_EPS)
    res = S21_MPI;
  else
    res = (S21_MPI / 2 - s21_asin(x));
  return res;
}

int check_nan(double x) { return (x != x); }

long double s21_asin(double x) {  // 3
  double tmp = x, result = x;
  if (x < -1 || x > 1) {
    result = S21_NAN;
  } else if (x == -1 || x == 1) {
    result = S21_MPI / 2 * x;
  } else {
    for (long double num = 1; s21_fabs(tmp) > S21_EPS; num++) {
      tmp *= ((x * x) * (2 * num - 1) * (2 * num - 1)) /
             ((2 * num) * (2 * num + 1));
      result += tmp;
    }
  }
  return result;
}

long double atan_range(double x) {
  long double result = x, temp = x, i = 1;
  while (s21_fabs(result) > S21_EPS) {
    result *= -1 * s21_pow(x, 2) * (2 * i - 1) / (2 * i + 1);
    i += 1;
    temp += result;
  }
  return temp;
}

long double s21_atan(double x) {  // 4
  long double temp = 0;
  /* (x больше 16 значащих цифр) или не определен NAN или INF */
  if ((x >= (double)LLONG_MAX || x <= (double)LLONG_MIN) || x != x) temp = x;
  /* x == 0,  atan равен 0 */
  else if (x == 0)
    temp = 0;
  /* x принадлежит (-1; 1), используется разложение в ряды Маклорена
  (частный случай, рядов Тейлора в окрестности нужной точки) */
  else if (x < 1 && x > -1)
    temp = atan_range(x);
  /* x == 1 , atan = PI / 4
  <https://www.rapidtables.com/math/trigonometry/arctan/arctan-of-1.html > */
  else if (x == 1)
    temp = S21_MPI / 4;
  /* x == -1, atan = -PI / 4 , т.к. atan нечетная функция поэтому
   * atan(-x)=-atan(x) */
  else if (x == -1)
    temp = -S21_MPI / 4;
  /* x > 0,  т.к. atan определен на промежутке от -PI/2 до PI/2 оно не может
   быть больше этих значений,
   чтобы максимально к ним приблизиться используется ряд Маклорена для точки 1/х
 */
  else if (x > 1)
    temp = S21_MPI / 2 - atan_range(1 / x);
  else if (x < -1)
    temp = -S21_MPI / 2 - atan_range(1 / x);
  return temp;
}

long double s21_ceil(double x) {  // 5
    long long int buf = (long long int)x;
    double result = (double)buf;
    /*
     * Если x - положительное целое, или отрицательное < -1, явное приведение типа в самом начале,
     * выполняет необходимое округление и мы, не заходя ни в одно из условий,
     * возвращаем результат этого приведения.
     */
    if (x >= (double)LLONG_MAX || x <= (double)LLONG_MIN || x != x) {  // Если x > 16 зн. или == inf/-inf/nan,
        result = x;  // возвращаем входящий аргумент.
    } else if (x > 0.0 && x != result) {  // Если x - положительное, с дробной частью,
        result += 1;  //  округляем до следующего целого.
    } else if (x < 0.0 && result == 0.0) {  // Если -1 < x < 0,
        result = -0.0;  // преобразуем 0.0 в -0.0.
    }
    return (long double)result;
}

long double s21_cos(double x) { return s21_sin(S21_MPI / 2 - x); }  // 6

long double s21_exp(double x) {  // 7
  long int i = 1;

  long double result = 1, degree = x, fact = 1;
  if (x > 709.7827) {
    result = S21_INFINITY;
  } else if (x < -708.3964) {
    result = 0.0;
  } else {
    while (i < 1000) {
      result += degree / fact;
      degree *= x;
      fact *= (i + 1);
      ++i;
    }
  }
  return result;
}

long double s21_fabs(double x) {  // 8
  if (x < 0.0) x = -x;
  return x;
}

long double s21_floor(double x) {  // 9
  double result = 0;
  if (x >= LLONG_MAX || x <= LLONG_MIN || x != x) {
    result = x;
  } else {
    long long int a = (long long int)x;
    double y = (double)a;
    if (x > 0.0 || x == y) {
      result = y;
    } else {
      result = y - 1;
    }
  }
  return result;
}

long double s21_fmod(double x, double y) {  // 10
  long double buf = S21_NAN;
  double del = 0.0;
  if (y != 0.0 || x != S21_INFINITY || x != -S21_INFINITY) {
    del = (long long int)(x / y);
    buf = x - del * y;
  }
  return buf;
}

// http://mat-an.ru/img/study/teylor.gif
long double s21_log(double x) {  // 11
  int i = 2;
  long double result = x - 1, degree = x - 1, tmp = x;
  if (x < 0) {
    result = S21_NAN;
  } else if (x == 0.0) {
    result = -S21_INFINITY;
  } else if (x < 2.0 + S21_EPS) {
    x = x - 1;

    while (s21_fabs(degree / i) > S21_EPS) {
      degree *= -1 * x;
      result += degree / i;
      ++i;
    }
  } else {  //  :\{ - ln(123,456) =ln(1,23456×10^2) = ln(1,23456)+ln(10^2) =
            //  = ln(1,23456)+2*ln(10) ≈ ln(1,23456)+2*2,3025851 \}
    int num = 0;
    while (tmp > 1.0) {
      tmp /= 10;
      num++;
    }
    result = s21_log(tmp) + num * S21_LN10;
  }
  return result;
}

long double s21_pow(double base, double exp) {  // 12
  long double number;
  if (base < 0) {
    if ((long int)exp == exp) {
      if (exp > 0) {
        number = base;
        for (long int i = 0; i < (long int)exp - 1; i++) number *= base;

      } else if (exp == 0) {
        number = 1;
      } else {
        number = 1 / base;
        for (long int i = 0; i < (long int)exp * (-1) - 1; i++) {
          number /= base;
        }
      }
    } else {
      if (exp == -S21_INFINITY || exp == S21_INFINITY) {
        if (base * (-1) < 1) {
          number = 0;
        } else if (base * (-1) == 1) {
          number = 1;
        } else {
          if (exp == -S21_INFINITY) {
            number = 0;
          } else {
            number = S21_INFINITY;
          }
        }
      } else {
        number = -S21_NAN;
      }
    }
  } else if (base == 0) {
    if (exp == 0)
      number = 1;
    else
      number = 0;
  } else if (base == 1) {
    number = 1;
  } else {
    if ((long int)exp == exp) {
      if (exp > 0) {
        number = base;
        for (long int i = 0; i < (long int)exp - 1; i++) {
          number *= base;
        }
      } else if (exp == 0) {
        number = 1;
      } else {
        number = 1 / base;
        for (long int i = 0; i < (long int)exp * (-1) - 1; i++) {
          number /= base;
        }
      }
    } else {
      number = s21_exp(exp * (double)s21_log(base));
    }
  }
  return number;
}

long double s21_sin(double x) {  // 13
  x = s21_fmod(x, 2 * S21_MPI);
  double n = x, sum = 0.0;
  int i = 1;

  do {
    sum += n;
    n *= -x * x / (2 * i * (2 * i + 1));
    i++;
  } while (s21_fabs(n) > S21_EPS);

  return sum;
}

long double s21_sqrt(double x) {  // 14
  long double result = 4;
  if (x < 0) {
    result = -S21_NAN;
  } else {
    long double temp = 0;
    while (s21_fabs(result - temp) > S21_EPS) {
      temp = result;
      result = (temp + x / temp) / 2;
    }
  }
  return result;
}

long double s21_tan(double x) { return s21_sin(x) / s21_cos(x); }  // 15
