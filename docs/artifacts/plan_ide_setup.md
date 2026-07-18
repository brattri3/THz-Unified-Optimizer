# План настройки IDE и виртуального окружения

В данном плане описываются шаги по инициализации виртуального окружения, установке необходимых зависимостей (включая `reportlab` и `lmfit`, используемые в скриптах проекта) и настройке конфигурации VS Code для корректной интеграции интерпретатора Python.

## Proposed Changes

### Окружение и зависимости

#### [MODIFY] [requirements.txt](file:///c:/WGP/requirements.txt)
Обновим файл зависимостей, добавив в него библиотеки, которые фактически используются в проекте (`reportlab` и `lmfit`), но отсутствовали в исходном `requirements.txt`.
Содержимое будет следующим:
```text
numpy>=1.20.0
scipy>=1.7.0
matplotlib>=3.4.0
lmfit>=1.2.0
reportlab>=4.0.0
```

#### [NEW] [settings.json](file:///c:/WGP/.vscode/settings.json)
Создадим файл конфигурации для VS Code, чтобы IDE автоматически подхватывала виртуальное окружение и правильно разрешала импорты из папки `unified_optimizer`.

Содержимое `.vscode/settings.json`:
```json
{
  "python.defaultInterpreterPath": "${workspaceFolder}/.venv/Scripts/python.exe",
  "python.terminal.activateEnvInTriggeredTerminal": true,
  "python.analysis.extraPaths": [
    "${workspaceFolder}/unified_optimizer"
  ]
}
```

### Выполнение команд (в процессе реализации после одобрения)
1. Создание виртуального окружения `.venv` на базе Python 3.14 через утилиту `uv`:
   ```bash
   uv venv --python 3.14
   ```
2. Установка зависимостей из обновленного `requirements.txt` во вновь созданное окружение:
   ```bash
   uv pip install -r requirements.txt
   ```

## Verification Plan

### Automated Tests
- Запуск проверки работоспособности окружения путем вызова интерпретатора из виртуального окружения и проверки импорта основных библиотек:
  ```bash
  .venv\Scripts\python.exe -c "import numpy, scipy, matplotlib, lmfit, reportlab; print('All imports successful!')"
  ```
- Прогонка существующих тестов, например:
  ```bash
  .venv\Scripts\python.exe -m unittest discover -s unified_optimizer -p "test_*.py"
  ```

### Manual Verification
- Проверка создания файлов конфигурации `.vscode/settings.json`.
- Запись изменений в `CHANGELOG.md` и отправка изменений в репозиторий через Git (`git add`, `git commit`, `git push`).
- Копирование артефактов в `docs/artifacts/`.
