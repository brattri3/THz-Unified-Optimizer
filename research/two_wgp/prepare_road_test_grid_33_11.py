#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""A11 — «дорога» файлов-заглушек для датасета test_grid_33_11 (P = 33 мкм, D = 11 мкм).

Что это
-------
Владелец измеряет на спектрометре, сохраняя трассы в журнальной нумерации
`NNN_aXX.txt` / `NNN_bg.txt`, и ведёт `NNN_META.txt` вручную. Этот скрипт
раскладывает заранее ПУСТЫЕ файлы в порядке запланированного прогона, чтобы во
время съёмки не приходилось помнить план: следующий незаполненный файл и есть
следующее действие. Ровно так была устроена дорога `data_pool/purewave/`.

Файлы нулевого размера `DataManager` пропускает штатно, поэтому недоснятая
дорога НЕ ломает анализ: обрабатывается ровно то, что успели снять.

План углов (обоснование — MEASUREMENT_ORDER.md, §2)
---------------------------------------------------
Порядок блоков выбран так, чтобы (а) при остановке в любой момент оставался
пригодный к анализу набор и (б) дрейф был отделим от физики:

  блок 1  яркий якорь 0° + тень B 80…100° через 2°   ← ядро измерения
  блок 2  яркая зона 10°-шагом, ±70°                  ← нормировка + гармоники
  блок 3  тень A −80…−100° через 2°                   ← КОНТРОЛЬНАЯ ГРУППА
  блок 4  тень B повторно, обратный ход (rep2)        ← чистая мера дрейфа
  блок 5  яркие якоря rep2 (+40, 0, −40)              ← контроль нормировки

Блоки 3 и 4 сняты рядом по времени ⇒ сравнение «тень A против тень B rep2» —
это π-периодичность при СОВПАДАЮЩЕМ времени (физика), а «тень B rep1 против
тень B rep2» — чистый дрейф. На purewave эти два эффекта оказались смешаны, и
разделить их удалось только постфактум.

Гардрейлы
---------
* создаются ТОЛЬКО пустые файлы; существующие файлы не трогаются НИКОГДА;
* если в целевом каталоге найден непустой файл с тем же именем — отказ (не
  перезапись), скрипт сообщает и выходит с кодом 2;
* по умолчанию dry-run, запись только с `--apply`.

Запуск:
    .venv\\Scripts\\python.exe research/two_wgp/prepare_road_test_grid_33_11.py
    .venv\\Scripts\\python.exe research/two_wgp/prepare_road_test_grid_33_11.py --apply
    .venv\\Scripts\\python.exe research/two_wgp/prepare_road_test_grid_33_11.py --status
"""
from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path

REPO = Path(__file__).resolve().parents[2]

# --------------------------------------------------------------- параметры образца
DATASET = "test_grid_33_11"     # имя зарезервировано в data_pool/DATASETS.md
P_UM = 33.0
D_UM = 11.0

# --------------------------------------------------------------- план углов
SHADOW_B = [80, 82, 84, 86, 88, 90, 92, 94, 96, 98, 100]
SHADOW_A = [-80, -82, -84, -86, -88, -90, -92, -94, -96, -98, -100]
BRIGHT = [10, 20, 30, 40, 50, 60, 70, -10, -20, -30, -40, -50, -60, -70]
ANCHORS_REP2 = [40, 0, -40]

BG_EVERY = 3                    # фон через каждые N сигнальных трасс


def build_plan():
    """-> список записей {seq, kind, angle, rep, note}. seq проставляется ниже."""
    items = []

    def sig(angle, rep, note):
        items.append({"kind": "sig", "angle": int(angle), "rep": rep, "note": note})

    # блок 1 — яркий якорь и тень B
    sig(0, 1, "блок 1: яркий якорь — нормировка до всего остального")
    for a in SHADOW_B:
        sig(a, 1, "блок 1: тень B (rep1), шаг 2°")
    # блок 2 — яркая зона
    for a in BRIGHT:
        sig(a, 1, "блок 2: яркая зона, шаг 10°")
    # блок 3 — тень A, контрольная группа для π-периодичности
    for a in SHADOW_A:
        sig(a, 1, "блок 3: тень A (rep1) — контроль π-периодичности")
    # блок 4 — тень B обратным ходом
    for a in reversed(SHADOW_B):
        sig(a, 2, "блок 4: тень B (rep2), обратный ход — мера дрейфа")
    # блок 5 — яркие якоря rep2
    for a in ANCHORS_REP2:
        sig(a, 2, "блок 5: якорь rep2 — контроль нормировки")

    # раскладка фонов и META по журнальной нумерации
    plan, seq, since_bg = [], 1, BG_EVERY
    plan.append({"seq": f"{seq:03d}", "kind": "meta", "angle": None, "rep": None,
                 "note": "старт сессии: заполнить шапку по шаблону 000_META_TEMPLATE.txt"})
    seq += 1
    for it in items:
        if since_bg >= BG_EVERY:
            plan.append({"seq": f"{seq:03d}", "kind": "bg", "angle": None, "rep": None,
                         "note": "фон (пустой пучок); покрывает следующие трассы до нового фона"})
            seq += 1
            since_bg = 0
        it["seq"] = f"{seq:03d}"
        plan.append(it)
        seq += 1
        since_bg += 1
    return plan


def fname(rec):
    if rec["kind"] == "meta":
        return f"{rec['seq']}_META.txt"
    if rec["kind"] == "bg":
        return f"{rec['seq']}_bg.txt"
    return f"{rec['seq']}_a{rec['angle']}.txt"


META_TEMPLATE = """\
# МЕТАДАННЫЕ (заполни через блокнот; это НЕ трейс измерения)
seq: 000
дата:
время (часы ПК):
температура зоны измерения, C:
влажность зоны, %:
PEC (кристалл 2-й гармоники) целевая T, C:
PEC реальная T в моменте, C:
# лазер (S401C/PM100A) — снимок мощности; рутинно ОДИН арм, 1-2 раза за сессию ОБА:
мощность перед генератором, мВт:
мощность перед детектором, мВт:
# --- ГЕОМЕТРИЯ УСТАНОВКИ ОБРАЗЦА (добавлено 2026-08-06 по урокам purewave) ---
# Без этих полей нельзя отличить «изменилась физика» от «съехал образец».
светлая апертура оправы, мм:
диаметр перетяжки пучка, мм:
образец стоит в перетяжке между линзами (да/нет):
положение пучка на апертуре (центр / смещение и куда, мм):
ЦЕНТРОВКА ПРОВЕРЕНА: пятно за полный поворот уходит не более чем на, мм:
образец переставлялся с прошлой сессии (да/нет; если да — что менялось):
окно сканирования, пс:
длительность одной трассы, мин:
комментарий:
"""


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--name", default=DATASET, help="имя датасета (каталог в data_pool/)")
    ap.add_argument("--apply", action="store_true", help="создать файлы (иначе dry-run)")
    ap.add_argument("--status", action="store_true", help="показать прогресс заполнения")
    args = ap.parse_args()

    out = REPO / "data_pool" / args.name
    plan = build_plan()
    n_sig = sum(1 for r in plan if r["kind"] == "sig")
    n_bg = sum(1 for r in plan if r["kind"] == "bg")

    if args.status:
        if not out.exists():
            print(f"каталог {out} не создан — сначала --apply")
            return 1
        done = todo = 0
        first_todo = None
        for rec in plan:
            f = out / fname(rec)
            if f.exists() and f.stat().st_size > 0:
                done += 1
            else:
                todo += 1
                if first_todo is None:
                    first_todo = rec
        print(f"{args.name}: заполнено {done} из {len(plan)}, осталось {todo}")
        if first_todo:
            print(f"следующий файл: {fname(first_todo)}  ({first_todo['note']})")
        return 0

    # --- проверка безопасности: ни одного непустого файла не трогаем
    clashes = []
    if out.exists():
        for rec in plan:
            f = out / fname(rec)
            if f.exists() and f.stat().st_size > 0:
                clashes.append(f.name)
    if clashes:
        print("ОТКАЗ: в каталоге уже есть НЕПУСТЫЕ файлы с плановыми именами "
              f"({len(clashes)} шт., например {clashes[:3]}).")
        print("Скрипт не перезаписывает данные. Выберите другое имя (--name) "
              "или разберитесь с каталогом вручную.")
        return 2

    print(f"датасет : {args.name}   (P = {P_UM} мкм, D = {D_UM} мкм, D/P = {D_UM / P_UM:.3f})")
    print(f"каталог : {out}")
    print(f"план    : {len(plan)} файлов = {n_sig} сигнальных + {n_bg} фонов + 1 META")
    print(f"время   : {(n_sig + n_bg) * 6.25 / 60:.1f} ч при 6:15/трасса (окно 130 пс)")
    print(f"          {(n_sig + n_bg) * 3.17 / 60:.1f} ч при 3:10/трасса (окно 65 пс) — рекомендуется")
    print()
    if not args.apply:
        for rec in plan[:8]:
            print(f"  {fname(rec):<18} {rec['note']}")
        print(f"  ... ещё {len(plan) - 8} файлов")
        print("\n(dry-run; повторите с --apply)")
        return 0

    out.mkdir(parents=True, exist_ok=True)
    created = skipped = 0
    for rec in plan:
        f = out / fname(rec)
        if f.exists():
            skipped += 1
            continue
        f.touch()
        created += 1
    (out / "000_META_TEMPLATE.txt").write_text(META_TEMPLATE, encoding="utf-8")

    plan_json = REPO / "research" / "two_wgp" / f"plan_{args.name}.json"
    plan_json.write_text(json.dumps(
        {"dataset": args.name, "P_um": P_UM, "D_um": D_UM,
         "D_over_P": round(D_UM / P_UM, 4), "bg_every": BG_EVERY,
         "n_signal": n_sig, "n_bg": n_bg, "plan": plan},
        ensure_ascii=False, indent=2), encoding="utf-8")

    print(f"создано {created} пустых файлов, пропущено {skipped} (уже были)")
    print(f"шаблон META : {out / '000_META_TEMPLATE.txt'}")
    print(f"машинный план: {plan_json}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
