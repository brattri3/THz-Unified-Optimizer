"""Переименование измерений аттенюатора под единую подпись `att-11-16` (санкция владельца 2026-08-04).

ЗАЧЕМ. Шесть наборов `356att`, `series1..5` — это ОДИН И ТОТ ЖЕ физический прибор
(«THz Tunable Precision Attenuator **ATT-11-16-CA85**», S/N #02721, апертура 85 мм; паспорт —
`attenuator_app/passports/ATT-11-16-CA85_02721.json`), снятый в разные дни и в разных условиях.
Имена `356att`/`seriesN` этого не сообщают и провоцируют трактовать наборы как разные образцы —
что уже влияло на выводы (A3 сравнивал их как «группу», не имея этого подтверждения в именах).

ЧТО НЕ ДЕЛАЕТСЯ И ПОЧЕМУ. Наборы НЕ сливаются в один датасет:
  * ключ хранения `DataManager` — тройка (датасет, угол, повтор) (`data_manager.py:47`), а
    `356att_0deg_rep1` и `series1_0deg_rep1` дали бы ОДИН ключ ⇒ второй файл молча затёр бы
    первый (`data_manager.py:53`) — тихая потеря данных;
  * это РАЗНЫЕ прямые модели: в `356att` вращается второй по лучу WGP (у детектора), в
    `series1-3` — первый; плёночные поляризаторы стоят по-разному (см. таблицу ниже).
Поэтому подпись прибора становится ПРЕФИКСОМ, а различитель прогона сохраняется.

КОНФИГУРАЦИИ (уточнено владельцем 2026-08-04, совпадает с `data_pool/two_wgp_attenuator/README.md`):
  * `356att` → `att-11-16-356`: плёночный поляризатор сразу после генератора, затем НЕподвижный
    WGP, затем вращающийся WGP; перед детектором плёночного поляризатора НЕТ.
  * `series1` → `att-11-16-s1`: плёночных поляризаторов НЕТ вовсе; вращается ПЕРВЫЙ WGP.
  * `series2` → `att-11-16-s2`, `series3` → `att-11-16-s3`: оба плёночных поляризатора на месте
    (вход и выход); вращается ПЕРВЫЙ WGP.
  * `series4`, `series5` → `att-11-16-s4-artifact`, `att-11-16-s5-artifact`: артефакты кругового
    алгоритма скана (по одному углу). Пометка вынесена В ИМЯ, чтобы их нельзя было случайно
    включить в фит — до сих пор это знание жило только в README.

БЕЗОПАСНОСТЬ. По умолчанию — ХОЛОСТОЙ прогон (ничего не меняется). Боевой прогон (`--apply`):
  1. отказывается работать, если рабочее дерево git грязное в `data_pool/` (нужна чистая точка отката);
  2. проверяет, что ни одно целевое имя не занято (коллизии целевых имён — отказ до начала);
  3. снимает SHA-256 каждого файла ДО, переименовывает через `git mv` (история файла сохраняется),
     снимает SHA-256 ПОСЛЕ и сверяет — содержимое обязано совпасть побайтно;
  4. сверяет через `DataManager`, что новый датасет отдаёт РОВНО те же углы и те же численные
     массивы, что старый отдавал до переименования (иначе — аварийный откат);
  5. пишет манифест `research/results/two_wgp/rename_manifest_<timestamp>.json`, по которому
     `--rollback <манифест>` возвращает всё обратно.

ВНИМАНИЕ — ЛОМАЮЩЕЕ ИЗМЕНЕНИЕ ВНЕ ЗОНЫ A. Имена датасетов зашиты в 43 файлах кода:
зона B (`research/experiments`, 20 файлов), продакшн (`scripts/`, `unified_optimizer/`, 11 —
ЗАМОРОЖЕН), зона A (6, чиню сам), зона C (`attenuator_app`, 3). Плюс 32 файла артефактов и
журналов ссылаются на старые имена — их править НЕЛЬЗЯ (append-only историческая запись),
для них служит таблица соответствия в `data_pool/two_wgp_attenuator/README.md`.
Поэтому `--apply` осмысленно запускать только вместе с правками зон B и C (хэндоффы) —
`--check-code` печатает точный список того, что сломается.

Запуск ИЗ КОРНЯ репозитория:
    .venv/Scripts/python.exe research/two_wgp/rename_attenuator_data.py            # холостой прогон
    .venv/Scripts/python.exe research/two_wgp/rename_attenuator_data.py --check-code
    .venv/Scripts/python.exe research/two_wgp/rename_attenuator_data.py --apply
    .venv/Scripts/python.exe research/two_wgp/rename_attenuator_data.py --rollback <манифест.json>
"""
from __future__ import annotations

import argparse
import hashlib
import json
import re
import subprocess
import sys
import time
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
REPO = HERE.parents[1]
sys.path.insert(0, str(REPO))

from unified_optimizer import config                              # noqa: E402
from unified_optimizer.data_manager import DataManager, parse_filename  # noqa: E402

DATA_DIR = REPO / "data_pool" / "two_wgp_attenuator"
OUT = REPO / "research" / "results" / "two_wgp"

DEVICE = "ATT-11-16-CA85 S/N #02721"
MAPPING = {
    "356att":  "att-11-16-356",
    "series1": "att-11-16-s1",
    "series2": "att-11-16-s2",
    "series3": "att-11-16-s3",
    "series4": "att-11-16-s4-artifact",
    "series5": "att-11-16-s5-artifact",
}
CONFIG_NOTE = {
    "att-11-16-356": "плёнка после генератора; WGP1 неподвижен, WGP2 (у детектора) вращается; "
                     "выходной плёнки НЕТ",
    "att-11-16-s1":  "плёнок НЕТ вовсе; вращается ПЕРВЫЙ WGP",
    "att-11-16-s2":  "обе плёнки (вход+выход); вращается ПЕРВЫЙ WGP",
    "att-11-16-s3":  "обе плёнки (вход+выход); вращается ПЕРВЫЙ WGP (повтор условий s2)",
    "att-11-16-s4-artifact": "АРТЕФАКТ кругового скана, 1 угол — не использовать в фитах",
    "att-11-16-s5-artifact": "АРТЕФАКТ кругового скана, 1 угол — не использовать в фитах",
}


# ------------------------------------------------------------------ утилиты
def sha256(path: Path) -> str:
    h = hashlib.sha256()
    h.update(path.read_bytes())
    return h.hexdigest()


def plan():
    """Список (старый путь, новый путь, старый датасет, новый датасет). Пустой = нечего делать."""
    rows = []
    for f in sorted(DATA_DIR.rglob("*.txt")):
        try:
            ds, angle, rep, kind = parse_filename(f.name)
        except ValueError:
            continue
        if ds not in MAPPING:
            continue
        m = re.match(r"^(.+)_([+-]?\d+deg_rep\d+_(?:sig|bg))$", f.stem)
        new_name = f"{MAPPING[ds]}_{m.group(2)}.txt"
        rows.append((f, f.with_name(new_name), ds, MAPPING[ds]))
    return rows


def git_clean(paths=("data_pool",)) -> bool:
    r = subprocess.run(["git", "status", "--porcelain", "--"] + list(paths),
                       cwd=REPO, capture_output=True, text=True)
    return r.returncode == 0 and not r.stdout.strip()


def snapshot(datasets):
    """Снимок данных глазами DataManager: углы + хеш числовых массивов на каждый угол."""
    dm = DataManager(config.DATA_DIR)
    snap = {}
    for ds in datasets:
        data = dm.get_data_for_dataset(ds)
        if not data:
            continue
        per = {}
        for angle, (t_s, E_s, t_b, E_b) in data.items():
            h = hashlib.sha256()
            for arr in (t_s, E_s, t_b, E_b):
                h.update(np.ascontiguousarray(arr, dtype=np.float64).tobytes())
            per[float(angle)] = h.hexdigest()
        snap[ds] = per
    return snap


# ------------------------------------------------------------- анализ кода
def check_code(verbose=True):
    """Где зашиты старые имена — что придётся починить вместе с переименованием."""
    pat = re.compile(r'["\'](' + "|".join(map(re.escape, MAPPING)) + r')["\']')
    zones = [("research/two_wgp", "A (моя — чиню сам)"),
             ("research/hypotheses", "A (моя — чиню сам)"),
             ("research/experiments", "B (ядро — ХЭНДОФФ)"),
             ("attenuator_app", "C (приложение — ХЭНДОФФ)"),
             ("docs/attenuator_app", "C (приложение — ХЭНДОФФ)"),
             ("unified_optimizer", "продакшн ЗАМОРОЖЕН — санкция владельца"),
             ("scripts", "продакшн ЗАМОРОЖЕН — санкция владельца"),
             ("research/paper", "D (статья — ХЭНДОФФ)")]
    files = subprocess.run(["git", "ls-files"], cwd=REPO,
                           capture_output=True, text=True).stdout.split("\n")
    agg = {}
    for rel in files:
        if not rel or Path(rel).suffix not in (".py", ".md", ".json"):
            continue
        if rel.startswith("research/results") or rel.startswith("research/logs"):
            zone = "артефакты/журналы — НЕ ПРАВИТЬ (история)"
        else:
            zone = "прочее"
            for pref, z in zones:
                if rel.startswith(pref):
                    zone = z
                    break
        try:
            n = len(pat.findall((REPO / rel).read_text(encoding="utf-8", errors="ignore")))
        except OSError:
            continue
        if n:
            a = agg.setdefault(zone, {"files": [], "hits": 0})
            a["files"].append((rel, n))
            a["hits"] += n
    if verbose:
        print(f"\nГде зашиты старые имена датасетов ({', '.join(MAPPING)}):")
        print(f"{'зона':<42}{'файлов':>8}{'вхожд.':>8}")
        for z, a in sorted(agg.items(), key=lambda x: -len(x[1]["files"])):
            print(f"{z:<42}{len(a['files']):>8}{a['hits']:>8}")
        for z, a in sorted(agg.items()):
            if z.startswith("артефакты"):
                continue
            print(f"\n  --- {z} ---")
            for rel, n in sorted(a["files"], key=lambda x: -x[1]):
                print(f"    {rel:<62}{n:>4}")
    return agg


# ------------------------------------------------------------ переименование
def apply(rows, verbose=True):
    if not git_clean():
        print("ОТКАЗ: рабочее дерево git грязное в data_pool/ — нужна чистая точка отката.")
        return None
    targets = [new for _, new, _, _ in rows]
    if len(set(targets)) != len(targets):
        print("ОТКАЗ: коллизия целевых имён (два файла претендуют на одно имя).")
        return None
    занятые = [t for t in targets if t.exists()]
    if занятые:
        print(f"ОТКАЗ: {len(занятые)} целевых имён уже существуют, напр. {занятые[0].name}")
        return None

    before_files = {str(old.relative_to(REPO)): sha256(old) for old, _, _, _ in rows}
    before_data = snapshot(sorted(MAPPING))

    moved = []
    for old, new, _, _ in rows:
        r = subprocess.run(["git", "mv", str(old.relative_to(REPO)), str(new.relative_to(REPO))],
                           cwd=REPO, capture_output=True, text=True)
        if r.returncode != 0:
            print(f"ОШИБКА git mv {old.name}: {r.stderr.strip()} — откатываю")
            for o, n in reversed(moved):
                subprocess.run(["git", "mv", str(n.relative_to(REPO)), str(o.relative_to(REPO))],
                               cwd=REPO, capture_output=True, text=True)
            return None
        moved.append((old, new))

    # 1. содержимое побайтно не изменилось
    bad = [new.name for old, new, _, _ in rows
           if sha256(new) != before_files[str(old.relative_to(REPO))]]
    # 2. DataManager отдаёт те же углы и те же массивы под новым именем
    after_data = snapshot(sorted(MAPPING.values()))
    mismatch = []
    for old_ds, new_ds in MAPPING.items():
        b, a = before_data.get(old_ds, {}), after_data.get(new_ds, {})
        if b != a:
            mismatch.append((old_ds, new_ds, len(b), len(a)))

    if bad or mismatch:
        print(f"ПРОВЕРКА НЕ ПРОЙДЕНА (хеши: {len(bad)}, данные: {len(mismatch)}) — откатываю")
        for o, n in reversed(moved):
            subprocess.run(["git", "mv", str(n.relative_to(REPO)), str(o.relative_to(REPO))],
                           cwd=REPO, capture_output=True, text=True)
        return None

    manifest = {
        "device": DEVICE, "when": time.strftime("%Y-%m-%d %H:%M:%S"),
        "session": "A", "nick": "A3", "mapping": MAPPING,
        "config_note": CONFIG_NOTE,
        "n_files": len(rows),
        "moves": [{"old": str(o.relative_to(REPO)), "new": str(n.relative_to(REPO)),
                   "sha256": sha256(n)} for o, n in moved],
        "verification": {
            "sha256_all_match": True,
            "datamanager_angles_and_arrays_identical": True,
            "datasets": {old: {"new": new, "n_angles": len(before_data.get(old, {}))}
                         for old, new in MAPPING.items()},
        },
    }
    OUT.mkdir(parents=True, exist_ok=True)
    path = OUT / f"rename_manifest_{time.strftime('%Y%m%d_%H%M%S')}.json"
    path.write_text(json.dumps(manifest, ensure_ascii=False, indent=2), encoding="utf-8")
    if verbose:
        print(f"\nПереименовано файлов: {len(moved)}. Проверки пройдены "
              f"(SHA-256 всех файлов совпали; DataManager отдаёт идентичные углы и массивы).")
        print(f"Манифест (для отката) -> {path}")
    return path


def rollback(manifest_path, verbose=True):
    m = json.loads(Path(manifest_path).read_text(encoding="utf-8"))
    ok = 0
    for mv in reversed(m["moves"]):
        new, old = REPO / mv["new"], REPO / mv["old"]
        if not new.exists():
            print(f"  пропуск: {mv['new']} не найден")
            continue
        r = subprocess.run(["git", "mv", mv["new"], mv["old"]], cwd=REPO,
                           capture_output=True, text=True)
        if r.returncode == 0:
            ok += 1
        else:
            print(f"  ОШИБКА: {mv['new']} -> {mv['old']}: {r.stderr.strip()}")
    if verbose:
        print(f"Откачено файлов: {ok} из {len(m['moves'])}")
    return ok


# -------------------------------------------------------------------- main
def main():
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0],
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--apply", action="store_true", help="боевой прогон (по умолчанию — холостой)")
    ap.add_argument("--check-code", action="store_true",
                    help="показать, где в коде зашиты старые имена")
    ap.add_argument("--rollback", metavar="МАНИФЕСТ", help="вернуть имена по манифесту")
    args = ap.parse_args()

    if args.rollback:
        sys.exit(0 if rollback(args.rollback) else 1)

    rows = plan()
    print(f"Прибор: {DEVICE}")
    print(f"Каталог: {DATA_DIR.relative_to(REPO)}")
    print(f"Файлов к переименованию: {len(rows)}\n")
    print(f"{'старый':<10}{'новый датасет':<26}{'файлов':>7}  конфигурация")
    for old_ds, new_ds in MAPPING.items():
        n = sum(1 for _, _, o, _ in rows if o == old_ds)
        print(f"{old_ds:<10}{new_ds:<26}{n:>7}  {CONFIG_NOTE[new_ds]}")
    print("\nПримеры (первые 4):")
    for old, new, _, _ in rows[:4]:
        print(f"    {old.name:<34} -> {new.name}")

    if args.check_code:
        check_code()

    if not args.apply:
        print("\nХОЛОСТОЙ прогон — ничего не изменено. Боевой запуск: --apply")
        print("ПОМНИТЬ: --apply ломает код, где имена зашиты (см. --check-code); "
              "запускать вместе с правками зон B и C.")
        return
    apply(rows)


if __name__ == "__main__":
    main()
