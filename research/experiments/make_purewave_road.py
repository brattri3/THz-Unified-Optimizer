"""Build a SEQUENCE-NUMBERED measurement road for PureWave so filenames sort in
the exact measuring order (owner request 2026-07-27), with interleaved METADATA
placeholders.

Simple names to type at the spectrometer:
  signal    ->  {seq:03d}_a{angle}.txt     e.g. 007_a82.txt, 015_a-40.txt, 002_a0.txt
  reference ->  {seq:03d}_bg.txt           e.g. 001_bg.txt
  METADATA  ->  {seq:03d}_META.txt         (pre-filled template; open in notepad, fill values)
Sorting by filename == the road order. META files carry date/time/temperature/
humidity/PEC-crystal-temperature snapshots for drift correlation; the parser
skips anything whose tag is not `a{int}` or `bg`, so META never breaks matching.

Road: passes over the angle set (alternating direction) so each angle's reps are
spread across the session; bg every ~3 signals; a META snapshot at start, every
~15 traces, and at the end.
"""
from pathlib import Path

POOL = Path(__file__).resolve().parents[2] / "data_pool" / "purewave"
POOL.mkdir(parents=True, exist_ok=True)

ZONES = {
    "neg_crossed": ([-95, -85], 1),
    "neg_cal":     ([-60, -40, -20], 2),
    "bright":      ([0, 20, 40, 55, 65, 72], 1),
    "crossed":     ([80, 82, 84, 86, 88, 90, 92, 94, 96, 98, 100], 3),
    "beyond":      ([105, 110], 2),
}
BG_EVERY = 3          # one reference per this many signals
META_EVERY = 15       # one metadata snapshot per this many traces (sig+bg)
MAX_REPS = max(r for _, r in ZONES.values())

META_TEMPLATE = (
    "# МЕТАДАННЫЕ (заполни через блокнот; это НЕ трейс измерения)\n"
    "seq: {seq:03d}\n"
    "дата: \n"
    "время (часы ПК): \n"
    "температура зоны измерения, C: \n"
    "влажность зоны, %: \n"
    "PEC (кристалл 2-й гармоники) целевая T, C: \n"
    "PEC реальная T в моменте, C: \n"
    "# лазер (S401C/PM100A) — снимок мощности; рутинно ОДИН арм, 1-2 раза за сессию ОБА:\n"
    "мощность перед генератором, мВт: \n"
    "мощность перед детектором, мВт: \n"
    "комментарий: \n"
)


def angles_for_pass(p):
    ang = sorted({a for (angles, reps) in ZONES.values() if reps >= p for a in angles})
    return ang if p % 2 == 1 else list(reversed(ang))


def main():
    removed = 0
    for f in POOL.glob("*.txt"):
        if f.stat().st_size == 0:
            f.unlink(); removed += 1
    # also clear previously-written META templates (they are non-empty but regenerated here)
    for f in POOL.glob("*_META.txt"):
        if f.exists():
            f.unlink(); removed += 1

    road = []
    seq = 0
    since_bg = 0
    since_meta = 0

    def add(kind, angle=None):
        nonlocal seq
        seq += 1
        if kind == "meta":
            tag = "META"; fn = f"{seq:03d}_{tag}.txt"
            (POOL / fn).write_text(META_TEMPLATE.format(seq=seq), encoding="utf-8")
        else:
            tag = "bg" if kind == "bg" else f"a{angle}"
            fn = f"{seq:03d}_{tag}.txt"
            (POOL / fn).write_text("", encoding="utf-8")
        road.append((seq, kind, angle, fn))

    def maybe_meta():
        nonlocal since_meta
        if since_meta >= META_EVERY:
            add("meta"); since_meta = 0

    add("meta")                                  # metadata snapshot at the very start
    add("bg")                                    # reference at the very start
    since_bg = 0; since_meta = 0
    for p in range(1, MAX_REPS + 1):
        for a in angles_for_pass(p):
            add("sig", a); since_bg += 1; since_meta += 1
            maybe_meta()
            if since_bg >= BG_EVERY:
                add("bg"); since_bg = 0; since_meta += 1
                maybe_meta()
    if since_bg > 0:
        add("bg")
    add("meta")                                  # closing metadata snapshot

    n_sig = sum(1 for r in road if r[1] == "sig")
    n_bg = sum(1 for r in road if r[1] == "bg")
    n_meta = sum(1 for r in road if r[1] == "meta")
    traces = n_sig + n_bg

    L = []
    L.append("# PureWave — дорога измерений (нумерованная, сортируется по имени)\n")
    L.append(f"**{n_sig} сигналов + {n_bg} фонов = {traces} трейсов**, плюс {n_meta} метаданных-заметок. "
             f"Один трейс 6:15 @130 пс.")
    L.append(f"- Время @6:15: **~{traces*6.25/60:.1f} ч**; при окне ~65 пс (×2): **~{traces*6.25/60*0.5:.1f} ч** "
             f"(тогда добавь повторы у скрещённого). Метаданные — ~1 мин каждая, в общий бюджет почти не идут.\n")
    L.append("## Имена файлов\n")
    L.append("- Сигнал:  `{seq}_a{угол}.txt`  (`007_a82.txt`, `015_a-40.txt`, `002_a0.txt`)")
    L.append("- Фон:     `{seq}_bg.txt`        (`001_bg.txt`)")
    L.append("- МЕТА:    `{seq}_META.txt`      (уже с шаблоном внутри — открой блокнотом и впиши значения)\n")
    L.append("Сортировка по имени = порядок. Про пары sig/bg и повторы не думай — свяжет скрипт потом.\n")
    L.append("## Что писать в META (шаблон уже внутри файла)\n")
    L.append("дата · время по часам ПК · температура зоны · влажность зоны · PEC (кристалл 2-й гармоники) "
             "целевая T и реальная T в моменте · комментарий. META стоит в начале, каждые ~15 трейсов и в конце "
             "— чтобы потом сопоставить дрейф сигнала с условиями.\n")
    L.append("## Дорога (по порядку)\n")
    L.append("| # | Что делать | Файл |")
    L.append("|---|---|---|")
    for seq_i, kind, angle, fn in road:
        if kind == "meta":
            what = "**МЕТА** — впиши условия (T/влажн./PEC/время)"
        elif kind == "bg":
            what = "ФОН (пустой пучок)"
        else:
            what = f"сигнал, угол {angle:+d}°"
        L.append(f"| {seq_i:03d} | {what} | `{fn}` |")
    L.append("")
    L.append("> Повторы каждого угла разнесены по проходам (чередование направления) — дрейф "
             "декоррелирует с углом. Скрещённое 80–100° — 3 повтора; негативы/beyond — 2; ярко — 1. "
             "Осталось время после короткого окна — продолжай нумерацию (лишние повторы скрещённого), "
             "скрипт подхватит их как rep4+.\n")

    L.append("## ЗАМЕТКА ДЛЯ ПОСТ-ОБРАБОТКИ (скрипт / другая сессия)\n")
    L.append("Разбор `{seq}_{tag}.txt`: tag=`a{целое}` → сигнал угла; tag=`bg` → референс; "
             "tag=`META` (или любой НЕ `a…`/`bg`) → метаданные, в трейс-обработке ПРОПУСТИТЬ.")
    L.append("Переименование в канон `purewave_{угол}deg_rep{N}_{sig,bg}.txt`:")
    L.append("1. Сортировать по seq. 2. Для сигнала: фон = ближайший `_bg` по seq (предпочт. предыдущий). "
             "3. `rep{N}` = номер появления этого угла среди сигналов (1-е → rep1). "
             "4. Сигнал → `_sig.txt`; скопировать фон → парный `_bg.txt` (дубли ок). "
             "5. Отрицательные углы — знаковое целое; DataManager понимает.")
    L.append("META-файлы собрать отдельно в таблицу условий (seq, дата, время, T зоны, влажн., PEC target/real) "
             "→ `research/results/validation/purewave_conditions.csv` для корреляции дрейфа. "
             "Готовый скрипт: `research/experiments/postproc_road.py` (сделать в отдельной сессии).\n")

    L.append("## ВАЖНО про геометрию\n")
    L.append("Паспортный период (PureWave P=25.5) ОТЛИЧАЕТСЯ от реального намотанного (погрешность намотки); "
             "истину даёт МИКРОФОТО. В GEOMETRY паспорт — приор, не ground truth.\n")

    (POOL / "MEASUREMENT_ORDER.md").write_text("\n".join(L), encoding="utf-8")

    print(f"removed {removed} old placeholders")
    print(f"road: {n_sig} sig + {n_bg} bg = {traces} traces + {n_meta} META notes")
    print(f"time ~{traces*6.25/60:.1f} h @6:15 | ~{traces*6.25/60*0.5:.1f} h @65ps")
    print("first 10:", [r[3] for r in road[:10]])


if __name__ == "__main__":
    main()
