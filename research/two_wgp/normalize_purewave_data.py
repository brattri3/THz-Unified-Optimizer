"""
Приведение имён трасс PureWave к конвенции DataManager.

КОНТЕКСТ
--------
Измерения PureWave лежат в `data_pool/purewave/` в «журнальной» нумерации прибора:
`NNN_aXX.txt` (сигнал под углом XX), `NNN_bg.txt` (фон), `NNN_META.txt` (метаданные).
`DataManager` их НЕ видит: он требует `{dataset}_{angle}deg_rep{N}_{sig|bg}.txt`.

Скрипт создаёт КАНОНИЧЕСКИЕ КОПИИ рядом с оригиналами — ровно как это уже сделано
для `data_pool/specac/` (там сырые `0_bg1`, `bg1`, ... сохранены, а копии названы
`specac_0deg_rep1_sig.txt`). Оригиналы НЕ переименовываются и НЕ удаляются:
журнальная нумерация — единственная связь трасс с файлами `NNN_META.txt`.

ДВА ДНЯ ИЗМЕРЕНИЙ (018_META.txt, комментарий владельца)
-------------------------------------------------------
- День 1 — 28.07.2026, seq 002..017: яркая зона, углы -95..80.
- День 2 — 03.08.2026, seq 019..043: скрещенная зона, углы 82..110 + обратный ход.
  «с этого файла начинается новый день измерений, 019_bg первое измерение в этот день,
   не нужно использовать этот файл для предыдущих измерений» ⇒ фон 019 НЕ распространяется
  назад на seq<=017. Это соблюдено в MAPPING: каждый угол берёт ПРЕДШЕСТВУЮЩИЙ фон
  своего дня.

ПОЛИТИКА ФОНА (data_pool/purewave/MEASUREMENT_ORDER.md)
-------------------------------------------------------
Фон снимается примерно каждые 2-3 угла и копируется в каждый слот, который покрывает.
Здесь фон приписан углам, ИЗМЕРЕННЫМ ПОСЛЕ НЕГО и до следующего фона.

РЕПЛИКИ
-------
Углы 94, 96, 98, 100, 105, 110 сняты дважды: на восходящем ходу (rep1) и на
нисходящем (rep2). Остальные углы — rep1.

Запуск (из корня репозитория):
    .venv\\Scripts\\python.exe research/two_wgp/normalize_purewave_data.py          # dry-run
    .venv\\Scripts\\python.exe research/two_wgp/normalize_purewave_data.py --apply
    .venv\\Scripts\\python.exe research/two_wgp/normalize_purewave_data.py --verify
"""
from __future__ import annotations

import argparse
import hashlib
import json
import sys
from datetime import datetime
from pathlib import Path

import numpy as np

REPO = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(REPO))

SRC_DIR = REPO / "data_pool" / "purewave"
DATASET = "purewave"
RESULTS = REPO / "research" / "results" / "two_wgp"

# (seq сигнала, угол, rep, seq фона, день)
MAPPING = [
    # --- день 1: 28.07.2026, META 001 -----------------------------------------
    ("003", -95, 1, "002", 1),
    ("004", -85, 1, "002", 1),
    ("005", -60, 1, "002", 1),
    ("007", -40, 1, "006", 1),
    ("008", -20, 1, "006", 1),
    ("009",   0, 1, "006", 1),
    ("011",  20, 1, "010", 1),
    ("012",  40, 1, "010", 1),
    ("013",  55, 1, "010", 1),
    ("015",  65, 1, "014", 1),
    ("016",  72, 1, "014", 1),
    ("017",  80, 1, "014", 1),
    # --- день 2: 03.08.2026, META 018 (граница дней) --------------------------
    ("020",  82, 1, "019", 2),
    ("021",  84, 1, "019", 2),
    ("022",  86, 1, "019", 2),
    ("024",  88, 1, "023", 2),
    ("025",  90, 1, "023", 2),
    ("026",  92, 1, "023", 2),
    ("028",  94, 1, "027", 2),
    ("029",  96, 1, "027", 2),
    ("030",  98, 1, "027", 2),
    ("032", 100, 1, "031", 2),
    ("033", 105, 1, "031", 2),
    ("035", 110, 1, "031", 2),   # META 034 (23:20) между 033 и 035, фона не было
    # обратный ход того же дня -> rep2
    ("037", 110, 2, "036", 2),
    ("038", 105, 2, "036", 2),
    ("039", 100, 2, "036", 2),
    ("041",  98, 2, "040", 2),
    ("042",  96, 2, "040", 2),
    ("043",  94, 2, "040", 2),
]

# исходные имена: сигнал `NNN_a<угол>.txt`, фон `NNN_bg.txt`
def src_sig(seq: str, angle: int) -> Path:
    return SRC_DIR / f"{seq}_a{angle}.txt"


def src_bg(seq: str) -> Path:
    return SRC_DIR / f"{seq}_bg.txt"


def dst(angle: int, rep: int, kind: str) -> Path:
    return SRC_DIR / f"{DATASET}_{angle}deg_rep{rep}_{kind}.txt"


def sha256(p: Path) -> str:
    h = hashlib.sha256()
    h.update(p.read_bytes())
    return h.hexdigest()


def check_sources() -> list[str]:
    """Проверки ДО любой записи: файлы на месте, непустые, сетка времени одна и та же."""
    problems: list[str] = []
    seen_dst: dict[Path, str] = {}
    t_ref = None
    for seq, angle, rep, bg_seq, _day in MAPPING:
        s, b = src_sig(seq, angle), src_bg(bg_seq)
        for p in (s, b):
            if not p.exists():
                problems.append(f"НЕТ ИСХОДНИКА: {p.name}")
                continue
            if p.stat().st_size == 0:
                problems.append(f"ПУСТОЙ ИСХОДНИК (заглушка): {p.name}")
                continue
            try:
                raw = np.loadtxt(p)
            except Exception as e:                       # noqa: BLE001
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
            if d in seen_dst and seen_dst[d] != srcp.name:
                problems.append(f"КОЛЛИЗИЯ ИМЁН: {d.name} <- {seen_dst[d]} и {srcp.name}")
            seen_dst[d] = srcp.name
    return problems


def plan() -> list[tuple[Path, Path]]:
    pairs: list[tuple[Path, Path]] = []
    for seq, angle, rep, bg_seq, _day in MAPPING:
        pairs.append((src_sig(seq, angle), dst(angle, rep, "sig")))
        pairs.append((src_bg(bg_seq), dst(angle, rep, "bg")))
    return pairs


def verify_datamanager() -> dict:
    """Загрузка датасета штатным DataManager + сверка массивов с исходниками."""
    from unified_optimizer import config
    from unified_optimizer.data_manager import DataManager

    dm = DataManager(config.DATA_DIR)
    keys = sorted(k for k in dm.raw_data if k[0] == DATASET)
    complete = [k for k in keys if "sig" in dm.raw_data[k] and "bg" in dm.raw_data[k]]

    mismatches = []
    for seq, angle, rep, bg_seq, _day in MAPPING:
        key = (DATASET, float(angle), rep)
        if key not in dm.raw_data:
            mismatches.append(f"DataManager не видит {key}")
            continue
        for kind, p in (("sig", src_sig(seq, angle)), ("bg", src_bg(bg_seq))):
            raw = np.loadtxt(p)
            E = raw[:, 1] - np.mean(raw[:50, 1])          # то же, что load_tds
            t_dm, E_dm = dm.raw_data[key][kind]
            if not (np.array_equal(t_dm, raw[:, 0]) and np.array_equal(E_dm, E)):
                mismatches.append(f"расхождение массива {key} {kind} <- {p.name}")

    data = dm.get_data_for_dataset(DATASET)
    return {
        "datasets_visible": dm.get_datasets(),
        "n_keys": len(keys),
        "n_complete_pairs": len(complete),
        "angles": sorted(data.keys()),
        "n_angles_unique": len(data),
        "reps_per_angle": {str(a): sum(1 for k in keys if k[1] == a) for a in sorted({k[1] for k in keys})},
        "array_mismatches": mismatches,
    }


def main() -> int:
    ap = argparse.ArgumentParser(description="Канонические имена для трасс PureWave")
    ap.add_argument("--apply", action="store_true", help="создать копии (по умолчанию dry-run)")
    ap.add_argument("--verify", action="store_true", help="только проверить уже созданные копии")
    args = ap.parse_args()

    print(f"Источник: {SRC_DIR}")
    print(f"В маппинге: {len(MAPPING)} сигнальных трасс, "
          f"{len({m[3] for m in MAPPING})} уникальных фонов, "
          f"{len({m[1] for m in MAPPING})} уникальных углов\n")

    if args.verify:
        info = verify_datamanager()
        print(json.dumps(info, ensure_ascii=False, indent=2))
        ok = (info["n_complete_pairs"] == len(MAPPING)) and not info["array_mismatches"]
        print("\n[OK] проверка пройдена" if ok else "\n[FAIL] см. выше")
        return 0 if ok else 1

    problems = check_sources()
    if problems:
        print("ПРОБЛЕМЫ С ИСХОДНИКАМИ (запись не выполнена):")
        for p in problems:
            print("  -", p)
        return 1
    print("[OK] исходники: все на месте, непустые, единая сетка времени, коллизий имён нет")

    pairs = plan()
    existing = [d for _s, d in pairs if d.exists()]
    if existing and not args.apply:
        print(f"\nВНИМАНИЕ: {len(existing)} целевых файлов уже существуют — будут перезаписаны.")

    if not args.apply:
        print("\nDRY-RUN. План (первые 8 из %d пар):" % len(pairs))
        for s, d in pairs[:8]:
            print(f"  {s.name:>16}  ->  {d.name}")
        print("  ...")
        print("\nЗапуск с --apply создаст копии.")
        return 0

    manifest = {
        "created": datetime.now().isoformat(timespec="seconds"),
        "dataset": DATASET,
        "mode": "copy (оригиналы сохранены, прецедент data_pool/specac)",
        "source_dir": str(SRC_DIR.relative_to(REPO)),
        "files": [],
    }
    for s, d in pairs:
        payload = s.read_bytes()
        d.write_bytes(payload)
        h_src, h_dst = sha256(s), sha256(d)
        assert h_src == h_dst, f"SHA-256 не совпал: {s.name} -> {d.name}"
        manifest["files"].append({"src": s.name, "dst": d.name, "sha256": h_dst})
    print(f"[OK] создано {len(pairs)} копий, SHA-256 совпал у всех")

    info = verify_datamanager()
    manifest["datamanager"] = info
    RESULTS.mkdir(parents=True, exist_ok=True)
    mpath = RESULTS / f"purewave_normalize_manifest_{datetime.now():%Y%m%d_%H%M%S}.json"
    mpath.write_text(json.dumps(manifest, ensure_ascii=False, indent=2), encoding="utf-8")

    print(f"[OK] DataManager: {info['n_complete_pairs']}/{len(MAPPING)} полных пар (sig+bg), "
          f"{info['n_angles_unique']} уникальных углов")
    if info["array_mismatches"]:
        print("[FAIL] расхождения массивов:")
        for m in info["array_mismatches"]:
            print("   -", m)
        return 1
    print("[OK] массивы совпали с исходниками бит-в-бит")
    print(f"Манифест: {mpath.relative_to(REPO)}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
