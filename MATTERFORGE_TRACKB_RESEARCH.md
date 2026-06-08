# MatterForge — двухдорожечный контур · RESEARCH-БЛОК Дорожки B (шаг 1)

> Статус: **research-блок завершён, ожидается подтверждение списка систем** (по порядку работ — шаг 2).
> Страницы НЕ генерируются до подтверждения. Это рабочие заметки шага 1, не контент страниц.

## Метод и честная оговорка по верификации
- Источники проверялись через `WebSearch` (сводки, построенные по самим страницам) + перекрёстная сверка
  между запросами. **`WebFetch` в этом окружении заблокирован (HTTP 403 на всех издателях, arXiv, PMC,
  Semantic Scholar)** — первоисточники-PDF напрямую не открывались. Где первоисточник на число
  засорсить не удалось — число НЕ ставится, помечено явно.
- Дифференциатор соблюдён строго: все кандидаты в B — **strongly correlated TM** (Fe/Cu/Mn-oxo центры,
  спиновые состояния), где **однодетерминантный DFT систематически врёт**. Не «общая квантовая химия»,
  не «в N раз точнее DFT».
- Что ИМЕННО валидируется (смысл плашки «ПОДТВЕРЖДЕНО ВАЛИДАЦИЕЙ МЕТОДА»):
  - для Fe/Cu-цеолитов — воспроизведение **экспериментально установленной электронной структуры /
    спин-состояния активного центра и его мультиреференсности** (качественно→полуколичественно);
  - **количественный** энергетический золотой стандарт сверки в этом режиме — **спин-состояния**
    (эксп. spin-crossover энтальпии / spin-forbidden переходы в ккал/моль).
  - Плашка НЕ означает «MatterForge изобрёл катализатор» и НЕ присваивает промышленные референции.

## Таблица-сводка (кандидаты Дорожки B)

| # | Система | Состав / тип / реакция | Эксп-данные для сверки | DOI (ключевой) | Active space | Кубиты (JW) | В B? |
|---|---------|------------------------|------------------------|----------------|--------------|-------------|------|
| 1 | **Fe-zeolite α-O** | высокоспиновый **S=2 Fe(IV)=O** (α-O) в Fe-ZSM-5/FER/BEA; α-Fe(II)+N₂O→α-O; **CH₄→CH₃OH** (RT), бензол→фенол (AlphOx, пилотировался) | эксп. основное состояние **S=2** Fe(IV)=O (MCD/Mössbauer/NRVS); мультиреференс oxyl↔ferryl (где функционалы DFT расходятся); DFT-барьер HAT из CH₄ ≈ 6.6 ккал/моль | 10.1038/nature19059 (Nature 2016); 10.1073/pnas.1721717115 (PNAS 2018); 10.1039/D5CS00496A (Chem Soc Rev 2026) | CAS(Fe 3d + O 2p) ≈ (8–12e, 8–9o) | **16–18** | **ДА** |
| 2 | **Cu-zeolite** | **[Cu₂O]²⁺** (Cu-ZSM-5) / **[Cu₃O₃]²⁺** (Cu-MOR), Cu(II)-oxo; **CH₄→CH₃OH** (низкоT, stepwise) | rRaman-полоса **22 700 см⁻¹** → бент моно-µ-оксо [Cu₂O]²⁺ (подтв. ¹⁸O-сдвигом, исключает прочие Cu-O); трёхъядерный [Cu₃O₃]²⁺ в MOR | 10.1073/pnas.0910461106 (PNAS 2009); 10.1038/ncomms8546 (Nat Commun 2015) | [Cu₂O]²⁺: 13o полн. / ~8–10o редуц. AVAS | **26 (полн.) / 16–20 (редуц.)**; [Cu₃O₃]²⁺ ~48 → нужен sub-cluster | **ДА*** |
| 3 | **Spin-state anchor** (метод-анкор, НЕ пром. витрина) | спин-состояния Fe(II/III/IV) и Mn/Co комплексов/порфиринов | **количественные эксп.** spin-crossover энтальпии / spin-forbidden переходы (ккал/моль); **DFT MUE >15 ккал/моль**, хим. точность 1 ккал/моль недостижима; CASPT2/NEVPT2/CCSD(T) — референс | 10.1039/c9cp00105k (PCCP 2019, Radoń); 10.3390/molecules28083487 (Molecules 2023, Por21); SSE17 (PMC11577268, DOI уточняется) | Fe 3d (+ double-shell) CAS(6e,5o)→(10e,12o) | **10 (мин., неточно) / 20–24 (точно)** | **ДА как анкор** |
| 4 | FeO⁺ газофаза | FeO⁺ + CH₄ → Fe⁺ + CH₃OH; парадигма two-state reactivity | эксп. **D₀(Fe⁺–O) НЕ засорсил чисто** (поиск путает с FeO**H**⁺ = 3.34–3.47 эВ); computed TS1 = 22.1 ккал/моль (расчёт, не эксп.) | 10.1021/ja971723u (JACS 1998); 10.1021/ja0017965 (JACS 2000) | Fe 3d + O 2p ≈ 8–9o | 16–18 | **КАНДИДАТ** (нет эксп. числа из первоисточника) |
| 5 | Cu/ZnO methanol | Cu/ZnO/Al₂O₃, синтез метанола; Cu(0)/Cu(I) **d¹⁰** | — | — | закрытая d¹⁰ → слабо коррелирован | — | **НЕТ** (вне strong-correlation) |
| 6 | Cr/SiO₂ Phillips | CrOₓ/SiO₂, полимеризация этилена; **ст. окисления Cr(II)/Cr(III) спорна** | нет единого чистого наблюдаемого (механизм инициации спорен) | — | — | — | **НЕТ / кандидат** (нет «известного ответа») |
| 7 | VPO | (VO)₂P₂O₇, n-бутан→малеиновый ангидрид; V(IV)/V(V) | DFT-механизмы есть, но нет чистого strong-correlation бенчмарка с эксп. подписью, отличного от DFT | 10.1021/ja3115746 (JACS 2013, DFT) | CAS(V 3d + O 2p) | ~16 | **НЕТ / кандидат** (доп. проработка) |

`*` — ДА с оговоркой: полные Cu-oxo кластеры выходят за коридор 12–20 кубитов; нужен редуцированный
AVAS / embedded sub-cluster (как уже делалось для NiO-18q и FeMoco).

## Полные ссылки
1. Snyder B.E.R., Vanelderen P., Bols M.L., Hallaert S.D., Böttger L.H., Ungur L., Pierloot K.,
   Schoonheydt R.A., Sels B.F., Solomon E.I. *The active site of low-temperature methane hydroxylation
   in iron-containing zeolites.* **Nature 536, 317–321 (2016).** DOI 10.1038/nature19059.
2. Snyder B.E.R. et al. *Structural characterization of a non-heme iron active site in zeolites that
   hydroxylates methane.* **PNAS (2018).** DOI 10.1073/pnas.1721717115.
3. *Alpha oxygen – a unique oxidation active site from a quantum chemical viewpoint.*
   **Chem. Soc. Rev. (2026).** DOI 10.1039/D5CS00496A.
4. Woertink J.S., Smeets P.J., Groothaert M.H., Vance M.A., Sels B.F., Schoonheydt R.A., Solomon E.I.
   *A [Cu₂O]²⁺ core in Cu-ZSM-5, the active site in the oxidation of methane to methanol.*
   **PNAS 106(45), 18908–18913 (2009).** DOI 10.1073/pnas.0910461106.
5. Grundner S. et al. *Single-site trinuclear copper oxygen clusters in mordenite for selective
   conversion of methane to methanol.* **Nat. Commun. 6, 7546 (2015).** DOI 10.1038/ncomms8546.
6. Radoń M. *Benchmarking quantum chemistry methods for spin-state energetics of iron complexes against
   quantitative experimental data.* **PCCP 21, 4854–4870 (2019).** DOI 10.1039/c9cp00105k.
7. Morgante P., Peverati R. *Comparison of the Performance of Density Functional Methods for the
   Description of Spin States and Binding Energies of Porphyrins (Por21).* **Molecules 28, 3487 (2023).**
   DOI 10.3390/molecules28083487.
8. *Performance of quantum chemistry methods for a benchmark set of spin-state energetics derived from
   experimental data of 17 transition metal complexes (SSE17).* PMC11577268. (DOI уточняется.)
9. Yoshizawa K. et al. *Methane–Methanol Conversion by MnO⁺, FeO⁺, and CoO⁺.* **JACS (1998).**
   DOI 10.1021/ja971723u. (FeO⁺ TS1 = 22.1 ккал/моль — расчёт.)
10. Shiota Y., Yoshizawa K. *Methane-to-Methanol Conversion by First-Row Transition-Metal Oxide Ions
    (ScO⁺…CuO⁺).* **JACS (2000).** DOI 10.1021/ja0017965.
11. *The Critical Role of Phosphate in VPO for n-Butane to Maleic Anhydride.* **JACS (2013).**
    DOI 10.1021/ja3115746.

## Рекомендация
- **Дорожка B (промышленные, валидировано):** №1 Fe-zeolite α-O + №2 Cu-zeolite.
  Обе — прямая функционализация метана (C–H активация) с **известным экспериментальным ответом** по
  электронной структуре активного центра; тематически усиливают флагман Track A (OCM, CH₄→C₂H₄):
  один движок, разный статус зрелости.
- **Метод-анкор (валидация движка, без промышленной витрины):** №3 спин-состояния
  (Radoń/SSE17/Por21) — единственная **количественная** эксп. сверка энергий в этом режиме.
- **Дисциплина (НЕ в B):** №5 Cu/ZnO (d¹⁰, не тот режим), №6 Cr-Phillips и №7 VPO
  (нет чистого «известного ответа»), №4 FeO⁺ (нет эксп. числа из первоисточника).

## Нужно подтверждение (шаг 2), затем — генерация
1. Финальный список систем Дорожки B (рекоменд.: №1, №2, +№3 как анкор).
2. Язык интерфейса страниц (RU / EN / билингва как у сайта).
3. Показывать ли блок «НЕ прошли в B» как прозрачность (рекоменд.: да).
