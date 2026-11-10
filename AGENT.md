# 1. Назначение проекта

Разработать поэтапную математическую и программную модель каталитической гидроизомеризации н-гексадекана на пористом цилиндрическом катализаторе Pt/SAPO-11.

Расчётное ядро реализовать на современном Fortran. Python использовать для:

- подготовки конфигураций;
- запуска расчётов;
- проверки входных данных;
- анализа результатов;
- визуализации;
- сравнения расчёта с экспериментом;
- идентификации параметров в будущих версиях.

Каждая версия V1–V6 обязана быть самостоятельно запускаемой и проверяемой.

Развитие:

**V1 — внутренняя реакционно-диффузионная модель гранулы**

→ **V2 — внешний массоперенос**

→ **V3 — гранула в реакторе**

→ **V4 — конечная цилиндрическая геометрия + тепловая модель**

→ **V5 — коксование и дезактивация**

→ **V6 — multiscale: гранула → кристалл SAPO-11 + итоговая валидация**

---

# 2. Общие требования

## 2.1. Языки

### Fortran

Использовать современный Fortran 2008/2018.

Fortran отвечает за:

- кинетику;
- расчёт коэффициентов переноса;
- пространственную дискретизацию;
- решение ODE/PDE;
- расчёт источников вещества и энергии;
- coupling между моделями.

Не реализовывать построение графиков внутри Fortran.

### Python

Python отвечает за:

- чтение пользовательского YAML/TOML-конфига;
- проверку единиц и обязательных параметров;
- генерацию входного файла Fortran;
- запуск executable;
- чтение результатов;
- графики;
- отчёты;
- автоматические сравнения;
- parameter sweep;
- sensitivity analysis;
- позднее — fitting параметров.

Основные библиотеки:

```text
numpy
pandas
matplotlib
scipy
pydantic
pyyaml
```

---

# 3. Принцип работы с параметрами

Каждый физический параметр должен иметь не только значение, но и происхождение.

Для параметра хранить:

```text
name
value
unit
source_type
source_reference
uncertainty
comment
```

`source_type`:

```text
measured
literature
estimated
fitted
temporary_demo
```

Пример:

```yaml
pellet:
  radius:
    value: 1.0e-3
    unit: m
    source_type: temporary_demo
    source_reference: ""
    uncertainty: null
```

## Режимы программы

### demo

Разрешает запуск с параметрами `literature`, `estimated`, `temporary_demo`.

В отчёте обязательно писать предупреждение:

```text
RESULTS ARE NOT EXPERIMENTALLY VALIDATED
```

### scientific

Запрещает использование `temporary_demo`.

Необходимые экспериментальные параметры должны быть получены или обоснованы литературой.

---

# 4. Единицы

Внутри Fortran использовать только SI.

Основные единицы:

```text
длина                 m
время                  s
температура            K
давление               Pa
концентрация           mol/m^3
молярный расход        mol/s
масса                  kg
плотность              kg/m^3
D                      m^2/s
теплопроводность       W/(m K)
Cp                     J/(kg K) либо J/(mol K), явно различать
энергия активации      J/mol
энтальпия реакции      J/mol
```

Нельзя смешивать минуты, kcal/mol, bar и SI внутри математического ядра.

Конвертация осуществляется на уровне ввода.

---

# 5. Химическая кинетика

## 5.1. Компоненты

Базовая модель содержит 15 компонент:

```text
C1  = n-C16H34
C2  = mono-iso-C16H34
C3  = multi-iso-C16H34
C4  = H2
C5  = C2H6-C5H12
C6  = C6H14-C14H30
C7  = n-C16H32
C8  = n-C16H33+
C9  = mono-iso-C16H33+
C10 = multi-iso-C16H33+
C11 = C2H4-C5H10
C12 = (C6H13-C14H29)+
C13 = mono-iso-C16H32
C14 = multi-iso-C16H32
C15 = C6H12-C14H28
```

В вычислительной реакционно-диффузионной модели основная переменная:

\[
C_i \quad [mol/m^3].
\]

Старую модель на мольных долях сохранить только как `legacy`-режим для воспроизведения статьи.

Кинетические константы, идентифицированные для мольных долей, запрещено автоматически считать эквивалентными константам концентрационной модели.

---

# 6. Схема реакций

Использовать текущую схему:

\[
1.\quad C_1 \rightleftharpoons C_7+C_4
\]

\[
2.\quad C_7+H^+ \rightleftharpoons C_8
\]

\[
3.\quad C_8\rightleftharpoons C_9
\]

\[
4.\quad C*8\rightleftharpoons C*{10}
\]

\[
5.\quad C*9\rightleftharpoons C*{10}
\]

\[
6.\quad C*8\rightarrow C*{11}+C\_{12}
\]

\[
7.\quad C*9\rightarrow C*{11}+C\_{12}
\]

\[
8.\quad C*{10}\rightarrow C*{11}+C\_{12}
\]

\[
9.\quad C*9\rightleftharpoons C*{13}+H^+
\]

\[
10.\quad C*{10}\rightleftharpoons C*{14}+H^+
\]

\[
11.\quad C*{12}\rightarrow C*{15}+H^+
\]

\[
12.\quad C\_{11}+C_4\rightarrow C_5
\]

\[
13.\quad C\_{13}+C_4\rightleftharpoons C_2
\]

\[
14.\quad C\_{14}+C_4\rightleftharpoons C_3
\]

\[
\boxed{
15.\quad C\_{15}+C_4\rightarrow C_6
}
\]

Стадия 15 является **необратимой**.

Обратимая стрелка, приведённая для этой стадии в исходной статье, является ошибкой, что подтверждено научным руководителем.

Следовательно:

- обратная стадия 15 в математической модели отсутствует;
- параметр \(k\_{-15}\) не используется;
- параметр \(E\_{-15}\) не используется;
- параметр \(k\_{0,-15}\) не используется;
- программный флаг для включения обратимости стадии 15 не требуется.

---

# 7. Скорости реакций

Текущая mass-action реализация:

\[
w*1=k_1C_1-k*{-1}C_7C_4
\]

\[
w*2=k_2C_7C*{H^+}-k\_{-2}C_8
\]

\[
w*3=k_3C_8-k*{-3}C_9
\]

\[
w*4=k_4C_8-k*{-4}C\_{10}
\]

\[
w*5=k_5C_9-k*{-5}C\_{10}
\]

\[
w_6=k_6C_8
\]

\[
w_7=k_7C_9
\]

\[
w*8=k_8C*{10}
\]

\[
w*9=k_9C_9-k*{-9}C*{13}C*{H^+}
\]

\[
w*{10}=k*{10}C*{10}-k*{-10}C*{14}C*{H^+}
\]

\[
w*{11}=k*{11}C\_{12}
\]

\[
w*{12}=k*{12}C\_{11}C_4
\]

\[
w*{13}=k*{13}C*{13}C_4-k*{-13}C_2
\]

\[
w*{14}=k*{14}C*{14}C_4-k*{-14}C_3
\]

\[
\boxed{
w*{15}=k*{15}C\_{15}C_4
}
\]

Стадия 15 всегда рассчитывается только в прямом направлении.

Концентрация \(H^+\) пока является модельным параметром.

Поддержать:

```text
acid_site_model = constant
acid_site_model = explicit       # future
```

Для `constant`:

\[
C*{H^+}=C*{H^+}^{const}.
\]

Не считать значение 1 физически доказанным.

---

# 8. Аррениусовская зависимость

Для каждой реакции:

\[
k*j(T)=k*{0,j}
\exp\left(-\frac{E_j}{RT}\right).
\]

Где:

\[
R=8.314462618\;J/(mol\,K).
\]

Для обратимых реакций:

\[
k*{-j}(T)=k*{0,-j}
\exp\left(-\frac{E\_{-j}}{RT}\right).
\]

Для стадии 15 рассчитывается только:

\[
\boxed{
k*{15}(T)=
k*{0,15}
\exp\left(-\frac{E\_{15}}{RT}\right)
}
\]

Обратная константа для стадии 15 отсутствует.

Для каждого \(k_0\) обязательно хранить его размерность.

Размерность зависит от порядка реакции.

---

# 9. Источники вещества

Для каждого компонента:

\[
R*i=\sum_j\nu*{ij}w_j.
\]

Для текущей модели:

\[
R_1=-w_1
\]

\[
R*2=w*{13}
\]

\[
R*3=w*{14}
\]

\[
R*4=w_1-w*{12}-w*{13}-w*{14}-w\_{15}
\]

\[
R*5=w*{12}
\]

\[
R*6=w*{15}
\]

\[
R_7=w_1-w_2
\]

\[
R_8=w_2-w_3-w_4-w_6
\]

\[
R_9=w_3-w_5-w_7-w_9
\]

\[
R*{10}=w_4+w_5-w_8-w*{10}
\]

\[
R*{11}=w_6+w_7+w_8-w*{12}
\]

\[
R*{12}=w_6+w_7+w_8-w*{11}
\]

\[
R*{13}=w_9-w*{13}
\]

\[
R*{14}=w*{10}-w\_{14}
\]

\[
R*{15}=w*{11}-w\_{15}.
\]

Код не должен жёстко дублировать эти уравнения во многих местах.

Хранить стехиометрическую матрицу:

\[
\nu\_{ij}
\]

и рассчитывать:

\[
R=\nu w.
\]

---

# 10. Численная архитектура

Использовать:

**Finite Volume Method + Method of Lines + stiff implicit solver.**

Основной backend:

```text
CVODE / BDF
```

Создать абстракцию:

```text
solver_api
```

чтобы позднее можно было использовать:

```text
CVODE
DLSODE
```

без изменения физической модели.

Параметры solver:

```text
rtol
atol
max_steps
max_step
initial_step
linear_solver
jacobian_mode
```

---

# V1. ВНУТРЕННЯЯ ДИФФУЗИЯ В ПОРИСТОЙ ЦИЛИНДРИЧЕСКОЙ ГРАНУЛЕ

# 11. Геометрия V1

Считать гранулу длинным цилиндром.

Использовать только радиальную координату:

\[
0\le r\le R_p.
\]

Предположение:

\[
L_p\gg 2R_p.
\]

Если оно не выполняется, V1 остаётся приближением до V4.

---

# 12. Основное уравнение V1

Для каждого компонента:

\[
\boxed{
\frac{\partial C*i}{\partial t}
=
\frac1r
\frac{\partial}{\partial r}
\left(
rD*{\mathrm{eff},i}
\frac{\partial C_i}{\partial r}
\right)

- R_i(C,T)
  }
  \]

При постоянном \(D\_{\mathrm{eff},i}\):

\[
\frac{\partial C*i}{\partial t}
=
D*{\mathrm{eff},i}
\left(
\frac{\partial^2C_i}{\partial r^2}

- \frac1r\frac{\partial C_i}{\partial r}
  \right)
- R_i.
  \]

  ***

# 13. Граничные условия V1

В центре:

\[
\boxed{
\left.
\frac{\partial C*i}{\partial r}
\right|*{r=0}=0
}
\]

На внешней поверхности:

\[
\boxed{
C*i(R_p,t)=C*{i,bulk}
}
\]

В V1 внешнее сопротивление отсутствует.

---

# 14. Начальные условия V1

Должны задаваться через конфиг:

\[
C*i(r,0)=C*{i,0}(r).
\]

Реализовать режимы:

```text
user_defined
bulk_equilibrated
hydrogen_pretreated
```

Не задавать одно начальное состояние как универсальное.

---

# 15. Пористость

Если известна apparent/particle density \(\rho*p\) и удельный объём пор \(V*{pore}\):

\[
\boxed{
\varepsilon*p=\rho_pV*{pore}
}
\]

при согласованных единицах.

Если известна skeletal density:

\[
\boxed{
\varepsilon*p=
\frac{V*{pore}}
{V\_{pore}+1/\rho_s}
}
\]

Не использовать bulk density всего каталитического слоя вместо плотности отдельной гранулы.

---

# 16. Характерный радиус пор

Для приближения цилиндрических пор:

\[
\boxed{
r*{pore}=
\frac{2V*{pore}}{S\_{BET}}
}
\]

где:

\[
V*{pore}=V*{micro}+V\_{meso}.
\]

Это эффективная оценка, а не реальное описание распределения пор.

---

# 17. Knudsen diffusion

Для газовой фазы:

\[
\boxed{
D*{K,i}
=
\frac{2}{3}r*{pore}
\sqrt{
\frac{8RT}{\pi M_i}
}
}
\]

где \(M_i\) в kg/mol.

---

# 18. Molecular diffusion

Архитектура должна поддерживать:

```text
diffusion_model = constant_deff
diffusion_model = knudsen
diffusion_model = bosanquet
diffusion_model = custom
```

Для смеси газов mixture-averaged приближение:

\[
D*{m,i}
=
\frac{1-y_i}
{\displaystyle\sum*{j\ne i}\frac{y*j}{D*{ij}}}.
\]

\(D\_{ij}\) должен вычисляться отдельным provider.

Можно реализовать Fuller correlation как одну из моделей, но она не должна быть единственным вариантом.

---

# 19. Bosanquet

При совместном влиянии молекулярной и Knudsen-диффузии:

\[
\boxed{
\frac1{D*{pore,i}}
=
\frac1{D*{m,i}}

- \frac1{D\_{K,i}}
  }
  \]

После поправки на пористость и извилистость:

\[
\boxed{
D*{\mathrm{eff},i}
=
\frac{\varepsilon_p}{\tau}
D*{pore,i}
}
\]

или только Knudsen:

\[
D*{\mathrm{eff},i}
=
\frac{\varepsilon_p}{\tau}D*{K,i}.
\]

Для жидкой или плотной фазы не применять Knudsen/Bosanquet автоматически.

Разрешить напрямую задавать измеренный:

\[
D\_{\mathrm{eff},i}.
\]

---

# 20. Tortuosity

\[
\tau
\]

является входным параметром.

Разрешить временную оценочную корреляцию только в `demo`.

Нельзя выдавать оценочную tortuosity за измеренную.

---

# 21. Модуль Тиле

Для тестовой реакции первого порядка:

\[
\boxed{
\phi=
R*p
\sqrt{
\frac{k}{D*{\mathrm{eff}}}
}
}
\]

Использовать только как диагностический параметр для задачи, для которой определение \(\phi\) корректно.

---

# 22. Аналитический effectiveness factor

Для длинного цилиндра и реакции первого порядка:

\[
\boxed{
\eta=
\frac{2I_1(\phi)}
{\phi I_0(\phi)}
}
\]

где \(I_0,I_1\) — модифицированные функции Бесселя.

Эту формулу использовать как:

1. диагностическую;
2. тест правильности численного решателя.

Не умножать реальные скорости PDE ещё раз на \(\eta\).

---

# 23. Численный effectiveness factor

Для произвольной нелинейной реакции \(j\):

\[
\boxed{
\eta*j=
\frac{
\displaystyle
\frac1{V_p}\int*{V*p}w_j(C)\,dV
}{
w_j(C*{surface})
}
}
\]

если знаменатель не равен нулю.

Для длинного цилиндра:

\[
\overline{R_i}
=
\frac{2}{R_p^2}
\int_0^{R_p}R_i(r)\,r\,dr.
\]

---

# 24. FVM V1

Разбить радиус на:

```text
N_r = 20 ... 200
```

контрольных объёмов.

На единицу длины:

\[
V*j=
\pi
\left(
r*{j+1/2}^2-r\_{j-1/2}^2
\right)
\]

\[
A*{j+1/2}=2\pi r*{j+1/2}.
\]

Поток через грань:

\[
N*{i,j+1/2}
=
-D*{\mathrm{eff},i}
\frac{C*{i,j+1}-C*{i,j}}
{\Delta r}.
\]

Баланс:

\[
\frac{dC*{i,j}}{dt}
=
-\frac{
A*{j+1/2}N\_{i,j+1/2}

- A*{j-1/2}N*{i,j-1/2}
  }{
  V_j
  }

* R\_{i,j}.
  \]

В центре площадь внутренней грани равна нулю.

---

# 25. Steady-state режим

Поддержать:

```text
simulation_mode = transient
simulation_mode = steady
```

Первую реализацию steady выполнить pseudo-transient способом:

интегрировать систему, пока:

\[
\max
\left|
\frac{dC*i}{dt}
\right|
<\varepsilon*{steady}.
\]

---

# 26. Что нужно получить для V1

## От химиков / лаборатории

Позднее запросить:

```text
радиус гранулы
длина гранулы
particle density
skeletal density, если имеется
BET surface area
micropore volume
mesopore volume
pore size distribution
tortuosity, если измерялась
температура
давление
состав исходной смеси
точная фаза реакционной среды
```

Главное:

```text
intrinsic kinetic parameters
```

полученные в условиях отсутствия существенных транспортных ограничений.

## Пока можно взять из литературы

```text
Rp
Lp
epsilon
tau
SBET
Vpore
rpore
Deff
```

но все они получают `source_type=literature` или `temporary_demo`.

---

# 27. Python-графики V1

Обязательно:

```text
C_i(r) для выбранных времен
C_i(t) в центре и у поверхности
average C_i(t)
reaction rates
η
φ для тестовой модели
conversion
grid-convergence plot
```

---

# 28. Проверка V1

Тест 1:

\[
R_i=0
\]

→ концентрация должна прийти к равномерному `bulk`.

Тест 2:

\[
D\_{\mathrm{eff}}\rightarrow\infty
\]

→ гранула становится практически однородной.

Тест 3:

простая первая реакция → численная \(\eta\) совпадает с аналитической цилиндрической формулой.

Тест 4:

сетки:

```text
Nr=25
Nr=50
Nr=100
Nr=200
```

должны показывать сходимость.

---

# V2. ВНЕШНИЙ МАССОПЕРЕНОС

# 29. Физический смысл

Теперь:

\[
C*{surface}\ne C*{bulk}
\]

в общем случае.

Между основным потоком и гранулой существует внешняя плёнка.

---

# 30. Новое граничное условие

На поверхности:

\[
\boxed{
D*{\mathrm{eff},i}
\left.
\frac{\partial C_i}{\partial r}
\right|*{R*p}
=
k*{m,i}
(C*{i,bulk}-C*{i,s})
}
\]

для выбранной ориентации координаты.

Реализация должна быть проверена по знаку потока отдельным unit test.

---

# 31. Sherwood

\[
Sh*i=
\frac{k*{m,i}d*{eq}}{D*{m,i}}.
\]

В качестве selectable correlation реализовать Wakao–Funazkri:

\[
\boxed{
Sh_i=
2+1.1Re_p^{0.6}Sc_i^{1/3}
}
\]

но correlation provider должен быть заменяемым.

---

# 32. Dimensionless numbers

\[
Re*p=
\frac{\rho_f u d*{eq}}{\mu}
\]

\[
Sc*i=
\frac{\mu}
{\rho_fD*{m,i}}.
\]

Тогда:

\[
\boxed{
k*{m,i}=
\frac{Sh_iD*{m,i}}{d\_{eq}}
}
\]

---

# 33. Эквивалентный диаметр

Для несферической гранулы поддержать:

\[
\boxed{
d\_{eq}=
\frac{6V_p}{A_p}
}
\]

и отдельный параметр:

\[
\psi
\]

— sphericity.

Correlation provider должен явно указывать, какой characteristic diameter использует.

---

# 34. Диагностика V2

Рассчитывать:

\[
C*{bulk}-C*{surface}
\]

и mass-transfer Biot-like number:

\[
Bi*{m,i}
=
\frac{k*{m,i}R*p}{D*{\mathrm{eff},i}}.
\]

---

# 35. Новые параметры V2

Нужно получить:

```text
скорость потока
плотность реакционной среды
динамическую вязкость
молекулярные коэффициенты диффузии
bed void fraction, если имеется слой
эквивалентный размер гранулы
sphericity
```

Позднее спросить:

- расход сырья;
- расход H2;
- давление;
- температуру;
- внутренний диаметр реактора;
- массу/объём слоя;
- способ упаковки катализатора.

---

# V3. СВЯЗЬ С FIXED-BED РЕАКТОРОМ

# 36. Новая модель

Теперь:

\[
C\_{bulk}
\]

не задаётся одинаковым во всём пространстве.

Он изменяется вдоль реактора.

Использовать координату:

\[
0\le z_R\le L_R.
\]

Не путать \(z_R\) с осевой координатой отдельной гранулы V4.

---

# 37. Базовая модель реактора

Первая V3:

```text
1D
steady
plug-flow
isothermal
```

Для компонента:

\[
\boxed{
\frac{dF_i}{dz_R}
=
A_R(1-\varepsilon_b)
\overline{R_i}
}
\]

где:

\[
F_i\quad [mol/s]
\]

\[
A_R\quad [m^2]
\]

\[
\varepsilon_b
\]

— межгранульная порозность слоя,

\[
\overline{R_i}
\]

— объёмно усреднённый источник внутри гранулы.

---

# 38. Bulk composition

\[
F_T=\sum_iF_i
\]

\[
y_i=\frac{F_i}{F_T}.
\]

Для ideal gas режима:

\[
\boxed{
C\_{i,bulk}
=
y_i\frac{P}{RT}
}
\]

Но thermodynamic provider должен быть заменяемым.

Не применять идеальный газ автоматически для жидкой/плотной фазы.

---

# 39. Coupling V3

На каждом \(z_R\):

```text
reactor state
    ↓
T, P, Cbulk
    ↓
pellet solver
    ↓
C(r)
    ↓
average Ri
    ↓
reactor balance
```

Реализовать сначала `pellet_mode=quasi_steady`.

То есть гранула для локального bulk-состава считается установившейся.

---

# 40. Pressure drop

Реализовать как переключаемую функцию.

Для packed bed базовая модель — Ergun:

\[
-\frac{dP}{dz_R}
=
\frac{
150\mu(1-\varepsilon_b)^2
}{
\varepsilon_b^3d_p^2
}u

- \frac{
  1.75\rho_f(1-\varepsilon_b)
  }{
  \varepsilon_b^3d_p
  }u^2.
  \]

Для цилиндрических частиц characteristic diameter и sphericity должны быть конфигурируемы.

---

# 41. WHSV

Поддержать ввод:

\[
\boxed{
WHSV=
\frac{\dot m*{feed}}{m*{cat}}
}
\]

но не приравнивать автоматически:

\[
t=1/WHSV.
\]

Это разные физические величины.

---

# 42. Показатели реактора

Обязательно:

\[
X*{C16}
=
1-
\frac{F_1}{F*{1,in}}.
\]

Селективность должна рассчитываться через отдельный observable provider.

Для lumped компонентов нельзя молча считать простую молярную селективность корректной.

Поддержать:

```text
molar
mass
carbon
paper_legacy
```

definitions.

---

# 43. Новые параметры V3

Спросить/найти:

```text
длина реактора
внутренний диаметр
масса катализатора
высота слоя
bed void fraction
feed molar/mass flow
H2/feed ratio
inlet pressure
inlet temperature
WHSV/LHSV
fluid density
viscosity
equation of state / phase information
```

Желательно получить экспериментальные профили или хотя бы:

```text
inlet state
outlet state
conversion
selectivity
```

для нескольких режимов.

---

# 44. Python V3

Графики:

```text
Fi(z)
yi(z)
conversion(z)
selectivity(z)
P(z)
surface vs bulk concentrations
pellet η versus reactor coordinate
```

---

# V4. КОНЕЧНЫЙ ЦИЛИНДР + ТЕПЛОВОЙ ЭФФЕКТ

# 45. Реальная геометрия гранулы

Перейти к осесимметричной 2D модели:

\[
C_i=C_i(r,z_p,t).
\]

Область:

\[
0\le r\le R_p
\]

\[
-\frac{L_p}{2}
\le z_p\le
\frac{L_p}{2}.
\]

Для экономии использовать симметрию:

\[
0\le z_p\le L_p/2.
\]

---

# 46. Массовое уравнение V4

\[
\boxed{
\frac{\partial C*i}{\partial t}
=
\frac1r
\frac{\partial}{\partial r}
\left(
rD*{\mathrm{eff},i}
\frac{\partial C_i}{\partial r}
\right)

- \frac{\partial}{\partial z*p}
  \left(
  D*{\mathrm{eff},i}
  \frac{\partial C_i}{\partial z_p}
  \right)
- R_i
  }
  \]

  ***

# 47. Граничные условия V4

Ось:

\[
\frac{\partial C_i}{\partial r}=0.
\]

Средняя плоскость:

\[
\frac{\partial C_i}{\partial z_p}=0.
\]

Боковая поверхность и торец:

\[
-D*{\mathrm{eff},i}\nabla C_i\cdot n
=
k*{m,i}(C*{i,s}-C*{i,b}).
\]

---

# 48. Уравнение энергии гранулы

\[
\boxed{
(\rho C*p)*{\mathrm{eff}}
\frac{\partial T*p}{\partial t}
=
\frac1r
\frac{\partial}{\partial r}
\left(
r\lambda*{\mathrm{eff}}
\frac{\partial T_p}{\partial r}
\right)

- \frac{\partial}{\partial z*p}
  \left(
  \lambda*{\mathrm{eff}}
  \frac{\partial T_p}{\partial z_p}
  \right)
- Q\_{\mathrm{rxn}}
  }
  \]

  ***

# 49. Источник тепла

\[
\boxed{
Q\_{\mathrm{rxn}}
=
-\sum_j\Delta H_j(T)w_j
}
\]

В единицах:

\[
W/m^3.
\]

---

# 50. Эффективная теплоёмкость

Поддержать прямой ввод:

\[
(\rho C*p)*{\mathrm{eff}}
\]

или оценочную модель:

\[
(\rho C*p)*{\mathrm{eff}}
=
(1-\varepsilon*p)\rho_sC*{p,s}

- \varepsilon*p\rho_fC*{p,f}.
  \]

  ***

# 51. Effective thermal conductivity

\[
\lambda\_{\mathrm{eff}}
\]

является отдельной transport property.

Реализовать:

```text
constant
correlation
custom
```

Не считать её автоматически равной теплопроводности чистого SAPO-11.

---

# 52. Тепловое граничное условие

\[
\boxed{
-\lambda\_{\mathrm{eff}}
\nabla T_p\cdot n
=
h(T_s-T_b)
}
\]

где \(h\) — внешний коэффициент теплоотдачи.

---

# 53. Nusselt

Поддержать correlation provider:

\[
Nu=
\frac{hd\_{eq}}{\lambda_f}.
\]

Как одну из моделей реализовать аналог packed-bed correlation:

\[
Nu=
2+1.1Re^{0.6}Pr^{1/3}
\]

с:

\[
Pr=
\frac{\mu C\_{p,f}}{\lambda_f}.
\]

---

# 54. Reactor energy coupling

Если V3 и V4 используются одновременно, тепло, покидающее гранулу, должно поступать в bulk-fluid energy balance.

Не учитывать одновременно одно и то же тепло:

1. через \(Q\_{rxn}\) в реакторном объёме;
2. и повторно через surface heat flux.

Должен использоваться один физически согласованный coupling.

---

# 55. Новые параметры V4

Нужно спросить/найти:

```text
точный радиус гранулы
точную длину гранулы
particle density
skeletal density
Cp catalyst
thermal conductivity catalyst
effective thermal conductivity pellet
Cp fluid
thermal conductivity fluid
fluid viscosity
fluid density
reaction enthalpies ΔHj
heat-transfer coefficient либо данные для его расчёта
reactor wall temperature
wall heat-transfer coefficient, если нужен
```

Особенно спросить химиков:

> Есть ли данные по тепловым эффектам стадий или следует рассчитывать их из термодинамических данных?

Для lumped/ionic компонентов не выдумывать \(\Delta H\).

---

# 56. Python V4

Графики:

```text
T(r,z)
Ci(r,z)
heat source map
temperature centre/surface
maximum ΔT pellet
reactor bulk temperature
heat flux
```

---

# V5. КОКСОВАНИЕ И ДЕЗАКТИВАЦИЯ

# 57. Главное требование

Не зафиксировать одну эмпирическую формулу коксования как истину.

Создать интерфейсы:

```text
coke_rate_model
activity_model
structure_degradation_model
```

---

# 58. Coke state

Ввести:

\[
q_c(r,z,t)
\]

например:

```text
kg_coke / kg_catalyst
```

или:

```text
mol_coke / m^3
```

Выбор должен быть явным.

---

# 59. Уравнение коксования

Общий интерфейс:

\[
\boxed{
\frac{\partial q*c}{\partial t}
=
r*{coke}(C,T,q_c;\theta_c)
}
\]

Для demo допускается модель вида:

\[
r\_{coke}
=
k_c(T)g(C)a^n
\]

\[
k*c=
k*{0,c}
\exp\left(
-\frac{E_c}{RT}
\right).
\]

Функция \(g(C)\) должна задаваться моделью и не выбираться Codex произвольно.

---

# 60. Activity

Поддержать несколько моделей, например:

\[
a=\exp(-\alpha q_c)
\]

или:

\[
a=
\frac1{1+\alpha q_c}.
\]

Но физическая модель выбирается только после появления данных.

Для реакции:

\[
\boxed{
w_j^{active}
=
a_jw_j^{intrinsic}
}
\]

Желательно поддержать разные активности:

```text
acid_activity
metal_activity
```

так как Pt и кислотные центры могут дезактивироваться по-разному.

---

# 61. Изменение структуры пор

Предусмотреть:

\[
\varepsilon_p=\varepsilon_p(q_c)
\]

\[
r*{pore}=r*{pore}(q_c)
\]

\[
\tau=\tau(q_c).
\]

После изменения структуры пересчитывать:

\[
D\_{\mathrm{eff},i}(q_c,T).
\]

В demo по умолчанию:

```text
structure_degradation = disabled
```

пока параметры отсутствуют.

---

# 62. Что нужно для V5

Запросить у химиков:

```text
conversion versus time-on-stream
selectivity versus time-on-stream
coke content versus operating time
TGA spent catalyst
BET fresh catalyst
BET spent catalyst
pore volume before/after
pore size distribution before/after
acid-site measurements before/after
Pt activity/dispersion before/after, если имеется
```

Нужно выяснить:

```text
какие компоненты считаются coke precursors
механизм дезактивации
обратима ли дезактивация
условия регенерации
```

Без этих данных V5 является архитектурной и демонстрационной моделью, а не прогнозной.

---

# 63. Python V5

Графики:

```text
coke(t)
activity(t)
conversion(t)
selectivity(t)
Deff(t)
porosity(t)
radial/2D coke distribution
fresh vs spent catalyst
```

---

# V6. MULTISCALE: ГРАНУЛА + КРИСТАЛЛ SAPO-11

# 64. Цель

Перестать представлять всю гранулу одной однородной пористой средой.

Разделить:

```text
межкристаллические мезо/макропоры
                ↓
кристалл SAPO-11
                ↓
микропоры
                ↓
активные центры
```

---

# 65. Pellet-scale field

Ввести:

\[
C_i^p(r,z,t).
\]

Диффузия:

\[
D\_{p,i}.
\]

Это транспорт по межкристаллической/мезопористой структуре гранулы.

---

# 66. Crystal-scale field

Для каждого representative crystal:

\[
C_i^c(\xi,t).
\]

Общее уравнение:

\[
\boxed{
\frac{\partial C*i^c}{\partial t}
=
\nabla*\xi
\cdot
\left(
D*{c,i}
\nabla*\xi C_i^c
\right)

- R_i^c
  }
  \]

где:

\[
D\_{c,i}
\]

— intracrystalline diffusivity.

Геометрия кристалла должна задаваться:

```text
slab
cylinder
sphere
custom
```

---

# 67. Межмасштабный обмен

Простейшее equilibrium boundary condition:

\[
\boxed{
C\_{i,surface}^c
=
K_iC_i^p
}
\]

где \(K_i\) — partition/adsorption coefficient.

Более полный вариант:

\[
N*i^{p\rightarrow c}
=
k*{int,i}
\left(
C*i^p-\frac{C*{i,s}^c}{K_i}
\right).
\]

---

# 68. Pellet-scale balance V6

\[
\boxed{
\frac{\partial C*i^p}{\partial t}
=
\nabla\cdot(D*{p,i}\nabla C_i^p)

- a\_{pc}N_i^{p\rightarrow c}

* R_i^p
  }
  \]

где:

\[
a\_{pc}
\]

— поверхность кристаллов на единицу объёма гранулы.

---

# 69. Расположение реакций

Каждая реакция должна иметь поле:

```text
reaction_location
```

например:

```text
acid_crystal
metal_surface
pellet
```

Это необходимо, чтобы не считать одну реакцию одновременно на двух масштабах.

Для Pt/SAPO-11 распределение металлических и кислотных центров должно задаваться из физической информации о катализаторе.

---

# 70. Что нужно получить для V6

Спросить/найти:

```text
размер кристаллов SAPO-11
распределение размеров
форма кристаллов
micropore dimensions
intracrystalline diffusivities
adsorption isotherms
partition coefficients
uptake curves
micro/mesopore volume fractions
intercrystalline porosity
Pt distribution
Pt particle size
location of Pt
acid-site concentration
acid-site distribution
```

Если доступны экспериментальные uptake/diffusion данные, использовать их вместо чисто оценочных корреляций.

---

# 71. Validation V6

Обязательно разделять:

```text
calibration dataset
validation dataset
```

Параметры запрещено подгонять по validation dataset.

Сравнивать:

```text
conversion
selectivity
product distributions
time dependence
temperature dependence
WHSV dependence
particle-size dependence
```

---

# 72. Метрики

Использовать:

\[
RMSE=
\sqrt{
\frac1N
\sum\_{k=1}^N
(y_k^{model}-y_k^{exp})^2
}.
\]

Также:

```text
MAE
relative error
maximum error
mass-balance residual
energy-balance residual
```

---

# 73. Parameter estimation

Python должен иметь отдельный модуль:

```text
fit/
```

который:

1. изменяет выбранные параметры;
2. запускает Fortran;
3. получает observables;
4. рассчитывает objective function;
5. использует scipy.optimize.

Нельзя смешивать fitting и обычный forward simulation.

---

# 74. Sensitivity analysis

Реализовать возможность анализировать:

\[
\frac{\partial Y}{\partial p_j}
\]

численно либо через finite differences.

Минимум:

```text
kinetic parameters
Deff
Rp
porosity
tortuosity
temperature
pressure
external km
coking parameters
```

---

# 75. Структура Fortran

Предлагаемая структура:

```text
src/
    kinds.f90
    constants.f90

    chemistry/
        species.f90
        stoichiometry.f90
        kinetics.f90
        arrhenius.f90

    transport/
        porosity.f90
        pore_geometry.f90
        molecular_diffusion.f90
        knudsen.f90
        effective_diffusion.f90
        external_mass_transfer.f90
        heat_transfer.f90

    pellet/
        pellet_1d.f90
        pellet_2d.f90
        pellet_sources.f90
        effectiveness.f90

    reactor/
        pfr.f90
        pressure_drop.f90
        reactor_energy.f90

    deactivation/
        coke.f90
        activity.f90

    multiscale/
        crystal.f90
        interscale_transfer.f90

    numerics/
        fvm_1d.f90
        fvm_2d.f90
        solver_api.f90
        cvode_backend.f90
        dlsode_backend.f90

    io/
        config.f90
        output.f90

    main.f90
```

---

# 76. Python

```text
python/
    config/
        models.py
        validation.py

    runner/
        build.py
        run.py

    analysis/
        observables.py
        balances.py
        convergence.py

    plots/
        pellet.py
        reactor.py
        temperature.py
        coke.py
        validation.py

    fit/
        objective.py
        optimize.py
        sensitivity.py
```

---

# 77. Выходные данные

Fortran должен сохранять:

```text
simulation metadata
parameter snapshot
solver statistics
time
coordinates
concentrations
reaction rates
Deff
surface fluxes
average sources
temperature
activity
coke
reactor states
mass-balance residual
energy-balance residual
```

Нельзя сохранять только итоговый conversion.

Нужна возможность воспроизвести расчёт полностью.

---

# 78. Обязательная проверка сохранения вещества

На каждом этапе рассчитывать elemental/material balance.

Проверять:

```text
carbon
hydrogen
```

настолько, насколько lumped reaction scheme позволяет это определить.

Если lumped species не позволяют однозначно восстановить атомный состав, выводить предупреждение.

---

# 79. Non-negativity

Концентрации физически:

\[
C_i\ge0.
\]

Нельзя просто после каждого шага делать:

```text
C=max(C,0)
```

и скрывать нестабильность.

Допускать только очень малые отрицательные значения порядка solver tolerance.

Существенные отрицательные значения должны приводить к diagnostic error.

---

# 80. Grid convergence

Для каждого пространственного решателя реализовать автоматический тест:

\[
E*N=
\frac{|Y*{2N}-Y*N|}
{\max(|Y*{2N}|,\epsilon)}.
\]

Расчёт считается пространственно разрешённым только после достижения заданного tolerance.

---

# 81. Solver validation

Для representative cases сравнить:

```text
CVODE BDF
DLSODE BDF
```

либо другой независимый stiff solver.

Результаты при уменьшении `rtol/atol` должны сходиться.

---

# 82. Regression tests между версиями

V2 при:

\[
k_m\rightarrow\infty
\]

должна воспроизводить V1.

V3 при постоянном bulk должна воспроизводить V2.

V4 при:

\[
L_p/R_p\rightarrow\infty
\]

и:

\[
T=\mathrm{const}
\]

должна воспроизводить V2/V3 pellet model.

V5 при:

\[
r\_{coke}=0
\]

должна воспроизводить V4.

V6 при:

\[
D\_{crystal}\rightarrow\infty
\]

и отсутствии interscale resistance должна приближаться к homogenized модели.

---

# 83. Ключевая научная оговорка

Имеющиеся сейчас кинетические параметры исходной статьи считаются **effective/apparent**, поскольку диффузионные эффекты явно не были отделены при их идентификации.

Поэтому запрещено интерпретировать расчёт:

```text
existing fitted kinetics
+
explicit diffusion
```

как окончательную физически валидированную модель.

Это допустимо только в `demo` для разработки программного комплекса.

Для итоговой `scientific` версии требуются intrinsic kinetic parameters.

---

# 84. Какой эксперимент нужен в будущем

Необходимо провести кинетические эксперименты с катализатором в условиях минимальных транспортных ограничений.

Желательно исследовать:

```text
несколько размеров частиц
несколько скоростей потока
несколько температур
малые/умеренные конверсии
```

Если дальнейшее:

```text
уменьшение размера частицы
```

не изменяет скорость реакции, внутреннее диффузионное ограничение мало.

Если дальнейшее:

```text
увеличение скорости внешнего потока
```

не изменяет скорость, внешнее сопротивление мало.

По этим данным требуется повторная идентификация:

\[
k\_{0,j},
\qquad
E_j.
\]

---

# 85. Таблица переходов между версиями

## V0 → V1

Нужно получить или оценить:

```text
Rp
Lp
porosity
pore volume
BET
pore radius
tortuosity
density
Mi
T
P
bulk composition
intrinsic kinetics
```

Результат:

```text
reaction + internal diffusion
```

## V1 → V2

Добавить:

```text
flow velocity
fluid density
fluid viscosity
Dm,i
particle equivalent diameter
sphericity
mass-transfer correlation
```

Результат:

```text
bulk → film → pellet
```

## V2 → V3

Добавить:

```text
reactor geometry
bed void fraction
catalyst mass
feed flow
H2/feed ratio
pressure
WHSV
pressure-drop data
```

Результат:

```text
fixed-bed reactor + pellet
```

## V3 → V4

Добавить:

```text
exact pellet radius
exact pellet length
Cp catalyst
lambda_eff
ΔH reaction
fluid Cp
fluid lambda
heat transfer
wall conditions
```

Результат:

```text
2D pellet + heat
```

## V4 → V5

Добавить:

```text
time-on-stream experiments
coke measurements
fresh/spent BET
fresh/spent pore distribution
activity decay
coke mechanism
```

Результат:

```text
dynamic catalyst degradation
```

## V5 → V6

Добавить:

```text
crystal size
crystal geometry
intracrystalline D
adsorption data
Pt distribution
acid-site distribution
micro/mesopore structure
independent validation experiments
```

Результат:

```text
multiscale validated digital twin
```

---

# 86. Что Codex НЕ должен делать

Codex запрещено:

1. придумывать отсутствующие физические параметры;
2. без предупреждения подменять intrinsic kinetic constants apparent-константами;
3. считать \(H^+=1\) физическим фактом;
4. вводить обратную реакцию для стадии 15;
5. вводить \(k*{-15}\), \(E*{-15}\) или \(k\_{0,-15}\);
6. использовать Knudsen diffusion без проверки выбранной phase model;
7. использовать bulk bed density как pellet density;
8. умножать PDE reaction rate дополнительно на analytical \(\eta\);
9. использовать slab effectiveness formula для цилиндра;
10. скрывать отрицательные концентрации clamp-операциями;
11. смешивать единицы;
12. считать fitting dataset одновременно validation dataset;
13. менять математическую модель для получения «красивых» результатов.

---

# 87. Порядок реализации Codex

Для каждой версии выполнять строго:

```text
1. математические структуры данных
2. входные параметры
3. физические функции
4. unit tests отдельных формул
5. пространственная дискретизация
6. solver
7. conservation tests
8. Python runner
9. графики
10. convergence tests
11. regression tests относительно предыдущей версии
12. документация
```

Не переходить к следующей версии, пока текущая не проходит тесты.

---

# 88. Итоговые критерии проекта

Финальная система должна позволять пользователю задать:

```text
структуру катализатора
размер гранулы
размер кристаллов
кинетические параметры
температуру
давление
расход
состав сырья
геометрию реактора
параметры массопереноса
тепловые параметры
параметры коксования
```

и получить:

```text
Ci(r,z,t)
Ci в кристалле
T(r,z,t)
Deff,i
reaction rates
η
surface fluxes
conversion
selectivity
cracking
coke
activity
reactor axial profiles
pressure drop
heat fluxes
mass/energy balance
```

с указанием происхождения всех использованных параметров и возможности полностью воспроизвести расчёт.

Конечная цель — не просто получить совпадающий график, а создать математически и физически прослеживаемую модель:

\[
\boxed{
\text{intrinsic kinetics}

- \text{transport}
- \text{geometry}
- \text{reactor}
- \text{heat}
- \text{deactivation}
- \text{multiscale structure}
  }
  \]

с независимой экспериментальной валидацией.
