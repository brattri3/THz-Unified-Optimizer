#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""A12 — нормализация имён датасета `test_grid_33_11` под конвенцию DataManager.

Прибор пишет журнальную нумерацию `NNN_aXX.txt` / `NNN_bg.txt`, `DataManager`
требует `test_grid_33_11_<angle>deg_rep<N>_<sig|bg>.txt`. Как в `purewave` и
`specac`, канонические имена появляются **копиями** рядом с оригиналами:
журнальная нумерация — единственная связь трасс с `NNN_META.txt`, поэтому
оригиналы не переименовываются и не удаляются.

Соответствие «журнальный номер → (угол, реплика)» НЕ придумывается здесь — оно
взято из плана, зафиксированного ДО измерения (`plan_test_grid_33_11.json`),
что исключает подгонку разметки под результат.

Правило фона: каждый `NNN_bg.txt` покрывает идущие ЗА ним сигнальные трассы до
следующего фона (та же политика, что заложена в маршрутный лист).

Запуск:
    .venv\\Scripts\\python.exe research/two_wgp/normalize_test_grid_33_11.py
    .venv\\Scripts\\python.exe research/two_wgp/normalize_test_grid_33_11.py --apply
    .venv\\Scripts\\python.exe research/two_wgp/normalize_test_grid_33_11.py --verify
"""
from __future__ import annotations

import argparse
import hashlib
import json
import shutil
import sys
from datetime import datetime
from pathlib import Path

import numpy as np

REPO = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(REPO))
DATASET = "test_grid_33_11"
SRC_DIR = REPO / "data_pool" / DATASET
PLAN_JSON = REPO / "research" / "two_wgp" / f"plan_{DATASET}.json"
RESULTS = REPO / "research" / "results" / "two_wgp"


def load_mapping():
    """-> [(seq, angle, rep, bg_seq)] по зафиксированному до измерения плану."""
    plan = json.loads(PLAN_JSON.read_text(encoding="utf-8"))["plan"]
    mapping, cur_bg = [], None
    for rec in plan:
        if rec["kind"] == "bg":
            cur_bg = rec["seq"]
        elif rec["kind"] == "sig":
            mapping.append((rec["seq"], int(rec["angle"]), int(rec["rep"]), cur_bg))
    return mapping


def src_sig(seq, angle):
    return SRC_DIR / f"{seq}_a{angle}.txt"


def src_bg(seq):
    return SRC_DIR / f"{seq}_bg.txt"


def dst(angle, rep, kind):
    return SRC_DIR / f"{DATASET}_{angle}deg_rep{rep}_{kind}.txt"


def sha256(p: Path) -> str:
    return hashlib.sha256(p.read_bytes()).hexdigest()


def check_sources(mapping):
    problems, seen, t_ref = [], {}, None
    for seq, angle, rep, bg_seq in mapping:
        if bg_seq is None:
            problems.append(f"НЕТ ФОНА перед сигналом {seq}_a{angle}")
            continue
        s, b = src_sig(seq, angle), src_bg(bg_seq)
        for p in (s, b):
            if not p.exists():
                problems.append(f"НЕТ ИСХОДНИКА: {p.name}")
                continue
            if p.stat().st_size == 0:
                problems.append(f"ПУСТОЙ (не снят): {p.name}")
                continue
            try:
                raw = np.loadtxt(p)
            except Exception as e:                                  # noqa: BLE001
                problems.append(f"НЕ ЧИТАЕТСЯ {p.name}: {e}")
                continue
            if raw.ndim != 2 or raw.shape[1] < 2:
                problems.append(f"ФОРМАТ не (время, поле): {p.name} shape={raw.shape}")
                continue
            if t_ref is None:
                t_ref = raw[:, 0]
            elif raw.shape[0] != len(t_ref) or not np.allclose(raw[:, 0], t_ref):
                problems.append(f"ДРУГАЯ СЕТКА ВРЕМЕНИ: {p.name}")
        for kind, srcp in (("sig", s), ("bg", b)):
            d = dst(angle, rep, kind)
            if d in seen and seen[d] != srcp.name:
                problems.append(f"КОЛЛИЗИЯ: {d.name} <- {seen[d]} и {srcp.name}")
            seen[d] = srcp.name
    return problems


def verify(mapping):
    from unified_optimizer import config
    from unified_optimizer.data_manager import DataManager

    dm = DataManager(config.DATA_DIR)
    dm.scan_directory()
    mism = []
    for seq, angle, rep, bg_seq in mapping:
        key = (DATASET, float(angle), rep)
        if key not in dm.raw_data:
            mism.append(f"DataManager не видит {key}")
            continue
        for kind, p in (("sig", src_sig(seq, angle)), ("bg", src_bg(bg_seq))):
            raw = np.loadtxt(p)
            E = raw[:, 1] - np.mean(raw[:50, 1])       # семантика load_tds
            t_dm, E_dm = dm.raw_data[key][kind]
            if not (np.array_equal(t_dm, raw[:, 0]) and np.array_equal(E_dm, E)):
                mism.append(f"расхождение массива {key} {kind} <- {p.name}")
    return {"datasets_visible": sorted(dm.get_datasets()),
            "keys": len([k for k in dm.raw_data if k[0] == DATASET]),
            "mismatches": mism}


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--apply", action="store_true")
    ap.add_argument("--verify", action="store_true")
    args = ap.parse_args()

    mapping = load_mapping()
    print(f"{DATASET}: {len(mapping)} сигнальных трасс по плану")

    if args.verify:
        res = verify(mapping)
        print(f"датасеты: {res['datasets_visible']}")
        print(f"ключей (угол, реплика): {res['keys']}")
        print("расхождений: " + (f"{len(res['mismatches'])}" if res["mismatches"] else "0 — OK"))
        for m in res["mismatches"][:10]:
            print("  " + m)
        return 1 if res["mismatches"] else 0

    problems = check_sources(mapping)
    if problems:
        print(f"ПРОБЛЕМЫ ({len(problems)}):")
        for p in problems[:20]:
            print("  " + p)
        if not args.apply:
            print("\n(dry-run; ошибки надо разобрать до --apply)")
        return 2

    pairs = []
    for seq, angle, rep, bg_seq in mapping:
        pairs.append((src_sig(seq, angle), dst(angle, rep, "sig")))
        pairs.append((src_bg(bg_seq), dst(angle, rep, "bg")))

    print(f"проверки исходников: OK; копий к созданию: {len(pairs)}")
    if not args.apply:
        for s, d in pairs[:4]:
            print(f"  {s.name} -> {d.name}")
        print(f"  ... ещё {len(pairs) - 4}")
        print("\n(dry-run; повторите с --apply)")
        return 0

    manifest = []
    for s, d in pairs:
        shutil.copy2(s, d)
        manifest.append({"src": s.name, "dst": d.name,
                         "sha256_src": sha256(s), "sha256_dst": sha256(d)})
    bad = [m for m in manifest if m["sha256_src"] != m["sha256_dst"]]
    RESULTS.mkdir(parents=True, exist_ok=True)
    stamp = datetime.now().strftime("%Y%m%d_%H%M%S")
    out = RESULTS / f"{DATASET}_normalize_manifest_{stamp}.json"
    out.write_text(json.dumps({"dataset": DATASET, "n_pairs": len(pairs),
                               "sha_mismatches": len(bad), "manifest": manifest},
                              ensure_ascii=False, indent=2), encoding="utf-8")
    print(f"создано копий: {len(pairs)}; несовпадений SHA-256: {len(bad)}")
    print(f"манифест: {out}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
