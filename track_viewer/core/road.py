# -*- coding: utf-8 -*-
"""Генерация «дороги» — заранее разложенных пустых файлов будущего прогона.

Зачем это нужно
---------------
Владелец измеряет, сохраняя трассы в журнальной нумерации, и во время съёмки не
должен помнить план: следующий незаполненный файл и есть следующее действие.
Файлы нулевого размера инструмент (и `DataManager` в ядре) пропускает штатно,
поэтому недоснятая дорога не ломает ни анализ, ни этот просмотрщик.

Обобщение, а не переписывание
-----------------------------
Основа — `research/two_wgp/prepare_road_test_grid_33_11.py` (зона A, читаем, не
правим). Оттуда взяты раскладка номеров, шаблон `META` и — главное — гардрейлы.
Захардкоженный план углов заменён параметрами, но так, чтобы значения по
умолчанию воспроизводили **ту же последовательность углов и фонов**, что лежит
в `data_pool/test_grid_33_11`. Это проверяется тестом, а не обещанием — иначе
«обобщение» незаметно разошлось бы с проверенной схемой.

Два отличия от того, что на диске, и оба намеренные:

* **нумерация сдвинута на единицу.** В старой схеме `001` занимал
  `001_META.txt`; по требованию владельца (ТЗ §2.1) `META` выведен из сквозной
  нумерации в единственный `META.txt` без номера, и первым файлом становится
  референс. Порядок углов от этого не меняется;
* **добавлен замыкающий фон.** ТЗ §5.1 требует «в начале, каждые n **и в
  конце**»; в старом скрипте последнего не было, из-за чего у последних трасс
  прогона нет позднего референса и режим «среднее двух» для них вырождается.

Почему блочный порядок — по умолчанию (решение владельца O-6)
-------------------------------------------------------------
Порядок обхода из `MEASUREMENT_ORDER.md` §2 решает две задачи сразу:

  блок 1  яркий якорь + тень (rep1)            ядро измерения
  блок 2  яркая зона крупным шагом             нормировка и гармоники
  блок 3  зеркальная тень                      контроль π-периодичности
  блок 4  та же тень обратным ходом (rep2)     чистая мера дрейфа
  блок 5  яркие якоря rep2                     контроль нормировки

Блоки 3 и 4 сняты рядом по времени, поэтому «зеркальная тень против тени rep2» —
это физика при совпадающем времени, а «тень rep1 против тени rep2» — чистый
дрейф. На `purewave` эти два эффекта оказались смешаны, и разделить их удалось
только постфактум. Монотонный обход по возрастанию такого разделения не даёт.

Гардрейлы (ТЗ §5.2) — не подлежат смягчению
-------------------------------------------
1. создаются ТОЛЬКО файлы нулевого размера;
2. существующий файл не перезаписывается НИКОГДА — ни пустой, ни заполненный;
3. найден непустой файл с плановым именем — операция отменяется ЦЕЛИКОМ;
4. сначала предпросмотр, создание — вторым явным действием;
5. удаление и усечение не выполняются ни при каких настройках.

Цена ошибки — затёртые измерения, которые нельзя восстановить.
"""
from __future__ import annotations

import io
import os

BG_START = "start"                  # один референс в начале
BG_START_END = "start_end"          # в начале и в конце
BG_EVERY_N = "every_n"              # в начале, каждые n трасс и в конце

ORDER_BLOCK = "block"
ORDER_ASCENDING = "ascending"

META_TEMPLATE = u"""\
# МЕТАДАННЫЕ (заполни через блокнот; это НЕ трейс измерения)
образец: {sample}
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


class Plan(object):
    """Параметры дороги. Значения по умолчанию = план `test_grid_33_11`."""

    __slots__ = ("sample", "coarse_from", "coarse_to", "coarse_step",
                 "fine_from", "fine_to", "fine_step", "fine_on",
                 "mirror_fine", "repeat_fine", "anchor", "anchors_rep2",
                 "bg_mode", "bg_every", "order", "repeats")

    def __init__(self, sample=u"образец",
                 coarse_from=-70, coarse_to=70, coarse_step=10,
                 fine_on=True, fine_from=80, fine_to=100, fine_step=2,
                 mirror_fine=True, repeat_fine=True,
                 anchor=0, anchors_rep2=(40, 0, -40),
                 bg_mode=BG_EVERY_N, bg_every=3,
                 order=ORDER_BLOCK, repeats=1):
        self.sample = sample
        self.coarse_from = int(coarse_from)
        self.coarse_to = int(coarse_to)
        self.coarse_step = int(coarse_step)
        self.fine_on = bool(fine_on)
        self.fine_from = int(fine_from)
        self.fine_to = int(fine_to)
        self.fine_step = int(fine_step)
        self.mirror_fine = bool(mirror_fine)
        self.repeat_fine = bool(repeat_fine)
        self.anchor = int(anchor)
        self.anchors_rep2 = tuple(int(a) for a in anchors_rep2)
        self.bg_mode = bg_mode
        self.bg_every = int(bg_every)
        self.order = order
        self.repeats = int(repeats)

    def validate(self):
        """-> список текстовых претензий; пустой список = план исполним."""
        bad = []
        if self.coarse_step <= 0:
            bad.append(u"шаг основного диапазона должен быть больше нуля")
        elif abs(self.coarse_to - self.coarse_from) < self.coarse_step:
            bad.append(u"основной диапазон уже одного шага")
        if self.fine_on:
            if self.fine_step <= 0:
                bad.append(u"шаг второго диапазона должен быть больше нуля")
            elif abs(self.fine_to - self.fine_from) < self.fine_step:
                bad.append(u"второй диапазон уже одного шага")
        if self.bg_mode == BG_EVERY_N and self.bg_every < 1:
            bad.append(u"референс не может ставиться реже чем через 1 трассу")
        if self.repeats < 1:
            bad.append(u"число повторов должно быть не меньше 1")
        if not str(self.sample).strip():
            bad.append(u"имя образца не задано — оно попадёт в META.txt и README.md")
        return bad


def _span(a, b, step):
    """Углы от a до b включительно с заданным шагом; направление — от a к b."""
    step = abs(int(step))
    out = []
    if a <= b:
        x = a
        while x <= b:
            out.append(int(x))
            x += step
    else:
        x = a
        while x >= b:
            out.append(int(x))
            x -= step
    return out


def build_plan(plan):
    """-> список записей {seq, kind, angle, rep, note}. `seq` — строка `NNN`."""
    bad = plan.validate()
    if bad:
        raise ValueError(u"; ".join(bad))

    coarse = _span(plan.coarse_from, plan.coarse_to, plan.coarse_step)
    fine = _span(plan.fine_from, plan.fine_to, plan.fine_step) if plan.fine_on else []
    mirror = [-a for a in fine] if (plan.fine_on and plan.mirror_fine) else []

    items = []

    def sig(angle, rep, note):
        items.append({"kind": "sig", "angle": int(angle), "rep": rep, "note": note})

    if plan.order == ORDER_ASCENDING:
        # Монотонный обход: меньше вращения оправы, но дрейф и угол смешаны.
        every = sorted(set(coarse + fine + mirror + [plan.anchor]))
        for rep in range(1, plan.repeats + 1):
            for a in every:
                sig(a, rep, u"обход по возрастанию%s"
                    % (u", повтор %d" % rep if plan.repeats > 1 else u""))
        return _lay_out(items, plan)

    # ------------------------------------------------------ блочный порядок
    # блок 1: яркий якорь — нормировка до всего остального, затем тень (rep1)
    sig(plan.anchor, 1, u"блок 1: яркий якорь — нормировка до всего остального")
    for a in fine:
        sig(a, 1, u"блок 1: тень (rep1), шаг %d°" % plan.fine_step)

    # блок 2: яркая зона. Положительные по возрастанию, затем отрицательные по
    # возрастанию модуля — так оправа возвращается к нулю один раз, а не дважды.
    rest = [a for a in coarse if a != plan.anchor]
    pos = sorted([a for a in rest if a > plan.anchor])
    neg = sorted([a for a in rest if a < plan.anchor], reverse=True)
    for a in pos + neg:
        sig(a, 1, u"блок 2: яркая зона, шаг %d°" % plan.coarse_step)

    # блок 3: зеркальная тень — контрольная группа для π-периодичности
    for a in mirror:
        sig(a, 1, u"блок 3: зеркальная тень (rep1) — контроль π-периодичности")

    # блок 4: та же тень обратным ходом — чистая мера дрейфа
    if plan.repeat_fine:
        for a in reversed(fine):
            sig(a, 2, u"блок 4: тень (rep2), обратный ход — мера дрейфа")

    # блок 5: яркие якоря rep2 — контроль нормировки
    for a in plan.anchors_rep2:
        sig(a, 2, u"блок 5: якорь rep2 — контроль нормировки")

    # Дополнительные повторы всего плана, если заказаны.
    if plan.repeats > 1:
        base = list(items)
        top = max([it["rep"] for it in base])
        for extra in range(2, plan.repeats + 1):
            for it in base:
                sig(it["angle"], it["rep"] + top * (extra - 1),
                    u"повтор %d всего плана" % extra)

    return _lay_out(items, plan)


def _lay_out(items, plan):
    """Расставляет референсы между сигнальными трассами и присваивает номера."""
    n = len(items)
    positions = set()
    if n:
        positions.add(0)                                    # референс в начале всегда
    if plan.bg_mode == BG_START_END and n:
        positions.add(n)
    elif plan.bg_mode == BG_EVERY_N:
        k = plan.bg_every
        i = k
        while i < n:
            positions.add(i)
            i += k
        positions.add(n)                                    # и в конце

    out = []
    seq = 1
    for i, it in enumerate(items):
        if i in positions:
            out.append({"seq": u"%03d" % seq, "kind": "bg", "angle": None, "rep": None,
                        "note": u"фон (пустой пучок); покрывает следующие трассы "
                                u"до нового фона"})
            seq += 1
        rec = dict(it)
        rec["seq"] = u"%03d" % seq
        out.append(rec)
        seq += 1
    if n in positions:
        out.append({"seq": u"%03d" % seq, "kind": "bg", "angle": None, "rep": None,
                    "note": u"замыкающий фон — контроль дрейфа за прогон"})
    return out


def fname(rec):
    if rec["kind"] == "bg":
        return u"%s_bg.txt" % rec["seq"]
    return u"%s_a%d.txt" % (rec["seq"], rec["angle"])


def summary(records):
    n_sig = len([r for r in records if r["kind"] == "sig"])
    n_bg = len([r for r in records if r["kind"] == "bg"])
    return {"n_sig": n_sig, "n_bg": n_bg, "n_files": n_sig + n_bg}


# ------------------------------------------------------------------ файлы
class Preview(object):
    """Что произойдёт при создании. Ничего не пишет на диск."""

    __slots__ = ("directory", "records", "to_create", "already_empty",
                 "conflicts", "extra_files")

    def __init__(self, directory, records):
        self.directory = directory
        self.records = records
        self.to_create = []      # имена, которых нет
        self.already_empty = []  # пустышки, которые останутся как есть
        self.conflicts = []      # НЕПУСТЫЕ файлы с плановыми именами — блокируют всё
        self.extra_files = []    # непустые файлы в каталоге вне плана

    @property
    def blocked(self):
        return bool(self.conflicts)

    def lines(self):
        out = [u"каталог: %s" % self.directory,
               u"план: %d файлов (%d трасс + %d фонов) + META.txt + README.md"
               % (len(self.records),
                  len([r for r in self.records if r["kind"] == "sig"]),
                  len([r for r in self.records if r["kind"] == "bg"]))]
        out.append(u"будет создано пустых файлов: %d" % len(self.to_create))
        if self.already_empty:
            out.append(u"уже лежат пустыми, останутся как есть: %d"
                       % len(self.already_empty))
        if self.extra_files:
            out.append(u"в каталоге есть посторонние непустые файлы (не трогаем): %d"
                       % len(self.extra_files))
        if self.conflicts:
            out.append(u"")
            out.append(u"ОТКАЗ: %d НЕПУСТЫХ файлов с плановыми именами."
                       % len(self.conflicts))
            out.append(u"Инструмент не перезаписывает измерения. Не создано ничего.")
            out.append(u"Выберите другой каталог или разберитесь с этим вручную.")
            for name in self.conflicts[:20]:
                out.append(u"    %s" % name)
            if len(self.conflicts) > 20:
                out.append(u"    … ещё %d" % (len(self.conflicts) - 20))
        return out


def preview(directory, records):
    """Проверка перед созданием. Диск только читается."""
    pv = Preview(directory, records)
    planned = set()
    if os.path.isdir(directory):
        existing = {}
        for name in os.listdir(directory):
            path = os.path.join(directory, name)
            if os.path.isfile(path):
                existing[name] = os.path.getsize(path)
        for rec in records:
            name = fname(rec)
            planned.add(name)
            if name not in existing:
                pv.to_create.append(name)
            elif existing[name] == 0:
                pv.already_empty.append(name)
            else:
                pv.conflicts.append(name)
        pv.extra_files = sorted([n for n, size in existing.items()
                                 if size > 0 and n not in planned])
    else:
        pv.to_create = [fname(r) for r in records]
    return pv


def apply(directory, records, plan, pv=None):
    """Создаёт пустые файлы, META.txt и README.md. -> (создано, пропущено).

    Пункты 2 и 3 гардрейлов ТЗ §5.2 исполняются здесь буквально: сначала полный
    предпросмотр, при единственном непустом конфликте — отказ до того, как
    создан хотя бы один файл; существующие файлы не открываются на запись
    вообще, поэтому усечь их невозможно даже по ошибке.
    """
    pv = pv or preview(directory, records)
    if pv.blocked:
        raise ValueError(u"в каталоге %d непустых файлов с плановыми именами; "
                         u"не создано ничего" % len(pv.conflicts))

    if not os.path.isdir(directory):
        os.makedirs(directory)

    created = skipped = 0
    for rec in records:
        path = os.path.join(directory, fname(rec))
        if os.path.exists(path):
            skipped += 1
            continue
        # 'xb' — создать или упасть. Не 'wb': тот усёк бы файл, появившийся
        # между предпросмотром и записью.
        try:
            with open(path, "xb"):
                pass
            created += 1
        except FileExistsError:
            skipped += 1

    _write_once(os.path.join(directory, "META.txt"),
                META_TEMPLATE.format(sample=plan.sample))
    _write_once(os.path.join(directory, "README.md"), readme_text(plan, records))
    return created, skipped


def _write_once(path, text):
    """Пишет файл, только если его ещё нет. Существующий не трогает никогда."""
    if os.path.exists(path):
        return False
    with io.open(path, "w", encoding="utf-8") as fh:
        fh.write(text)
    return True


def readme_text(plan, records):
    s = summary(records)
    lines = [
        u"# %s — угловой прогон" % plan.sample,
        u"",
        u"**Статус: дорога разложена, измерения не начаты.** Все трассы — файлы нулевого",
        u"размера; инструменты их пропускают, поэтому недоснятый прогон обрабатывается",
        u"штатно: считается ровно то, что успели снять.",
        u"",
        u"Разложено `track_viewer` (генератор трека).",
        u"",
        u"## План",
        u"",
        u"| Параметр | Значение |",
        u"|---|---|",
        u"| основной диапазон | %+d…%+d°, шаг %d° |" % (
            plan.coarse_from, plan.coarse_to, plan.coarse_step),
    ]
    if plan.fine_on:
        lines.append(u"| второй диапазон | %+d…%+d°, шаг %d° |" % (
            plan.fine_from, plan.fine_to, plan.fine_step))
        lines.append(u"| зеркальный второй диапазон | %s |"
                     % (u"да" if plan.mirror_fine else u"нет"))
        lines.append(u"| повтор второго диапазона обратным ходом | %s |"
                     % (u"да" if plan.repeat_fine else u"нет"))
    lines += [
        u"| яркий якорь | %+d° |" % plan.anchor,
        u"| порядок обхода | %s |" % (u"блочный" if plan.order == ORDER_BLOCK
                                      else u"по возрастанию"),
        u"| референсы | %s |" % {
            BG_START: u"один в начале",
            BG_START_END: u"в начале и в конце",
            BG_EVERY_N: u"в начале, каждые %d трасс и в конце" % plan.bg_every,
        }[plan.bg_mode],
        u"| всего файлов | %d = %d трасс + %d фонов |" % (
            s["n_files"], s["n_sig"], s["n_bg"]),
        u"",
        u"## Порядок съёмки",
        u"",
        u"Следующий незаполненный файл и есть следующее действие.",
        u"",
        u"| Файл | Что снимать |",
        u"|---|---|",
    ]
    for rec in records:
        what = (u"фон (пустой пучок)" if rec["kind"] == "bg"
                else u"угол %+d°%s" % (rec["angle"],
                                       u", повтор %d" % rec["rep"]
                                       if rec["rep"] and rec["rep"] > 1 else u""))
        lines.append(u"| `%s` | %s |" % (fname(rec), what))
    lines += [
        u"",
        u"## Перед стартом",
        u"",
        u"Заполнить `META.txt` — без полей геометрии установки нельзя отличить",
        u"«изменилась физика» от «съехал образец». Окно сканирования выбрать один раз",
        u"и не менять в течение прогона.",
        u"",
    ]
    return u"\n".join(lines)
