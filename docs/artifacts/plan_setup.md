# План исправления импортов и путей для запуска THz Unified Optimizer

Этот план направлен на устранение ошибок импорта (`ModuleNotFoundError: No module named 'unified_optimizer'`) и неверных путей к данным, чтобы сделать возможным ручной запуск фреймворка и его тестов.

## User Review Required

> [!IMPORTANT]
> Мы переименуем папку исходного кода `src` в `unified_optimizer`, чтобы импорты вида `from unified_optimizer...` работали "из коробки" при запуске из корня проекта. Также будут обновлены пути к данным в `config.py` и инструкции в `README.md`.

## Proposed Changes

### Реструктуризация проекта и конфигурация

---

#### [NEW] [unified_optimizer](file:///c:/THz-Unified-Optimizer/unified_optimizer)
Переименование директории `src` в `unified_optimizer`.

#### [DELETE] [src](file:///c:/THz-Unified-Optimizer/src)
Удаление старого имени директории.

#### [MODIFY] [config.py](file:///c:/THz-Unified-Optimizer/src/config.py)
Исправление путей к данным и результатам:
- Изменение `DATA_DIR = BASE_DIR / "data_pool"` на `DATA_DIR = BASE_DIR.parent / "data_pool"`.
- Изменение `RESULTS_DIR = BASE_DIR / "results"` на `RESULTS_DIR = BASE_DIR.parent / "results"` (или сохранение внутри пакета, но лучше вынести результаты в корень проекта для удобства).

#### [MODIFY] [README.md](file:///c:/THz-Unified-Optimizer/README.md)
Обновление путей структуры проекта и инструкций по запуску с учетом нового имени директории `unified_optimizer`.

## Verification Plan

### Automated Tests
- Запуск тестов из виртуального окружения:
  ```powershell
  .venv\Scripts\python -m unittest discover -s unified_optimizer/tests
  ```

### Manual Verification
- Запуск основного пайплайна для проверки работы оптимизатора и генерации отчетов:
  ```powershell
  $env:PYTHONPATH="."
  .venv\Scripts\python unified_optimizer/run_pipeline.py
  ```
- Проверка создания отчета в папке результатов.
