# Hydroisomerization

Поэтапная математическая и программная модель гидроизомеризации
н-гексадекана на пористом цилиндрическом катализаторе Pt/SAPO-11.

## V1: внутренняя диффузия в грануле

Текущая реализация включает:

- современный Fortran и 15-компонентную кинетическую схему;
- стехиометрическую матрицу и необратимую стадию R15;
- Arrhenius-параметры статьи с переводом энергии в SI;
- Knudsen-оценку эффективной диффузии;
- радиальный finite volume method для длинного цилиндра;
- Method of Lines и CVODE/BDF;
- Python-валидацию YAML, запуск, уникальные каталоги результатов и графики;
- контроль отрицательных концентраций без clamp;
- FVM-балансовую невязку и тест сходимости сетки.

Текущая конфигурация является демонстрационной:

```text
RESULTS ARE NOT EXPERIMENTALLY VALIDATED
```

Кинетика статьи является apparent и идентифицирована по мольным долям и
условному времени контакта. Её соединение с явной диффузией не является
окончательной физически валидированной моделью.

## Сборка

Требуется MSYS2 UCRT64 с GNU Fortran, CMake, Ninja и SUNDIALS.

```powershell
cmake -S . -B build -G Ninja `
  -DCMAKE_BUILD_TYPE=Debug `
  -DCMAKE_PREFIX_PATH=C:/msys64/ucrt64
cmake --build build
ctest --test-dir build --output-on-failure
```

## Запуск

```powershell
python python/run_v1.py config.yaml
```

Изменение одного или нескольких параметров без редактирования YAML:

```powershell
python python/run_v1.py config.yaml `
  --set conditions.temperature.value=613.15 `
  --set pellet.radius.value=1.0e-3
```

Каждая уникальная конфигурация и версия исходного кода получает отдельную
папку в `results/`. Повторный запуск не перезаписывает предыдущие результаты.

Проверка сеточной сходимости:

```powershell
python python/run_v1.py config.yaml --grid-convergence
```

Python-тесты:

```powershell
python -m unittest discover -s python/tests -v
```

## Результаты

Для каждого запуска сохраняются:

- нормализованный конфиг и Fortran namelist;
- радиальные концентрационные профили;
- объёмно-усреднённые концентрации, скорости и источники;
- поверхностные потоки;
- effectiveness factors;
- FVM-балансовые невязки;
- статистика CVODE;
- conversion и PNG-графики;
- hash конфигурации и исходного кода.

Файлы `kinetic_model.py` и `kinetic_model_advanced.py` сохранены как legacy-
воспроизведение исходной модели статьи.
