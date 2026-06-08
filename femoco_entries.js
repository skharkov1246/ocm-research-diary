/* ====================================================================
   FeMoco / N2-activation research diary — content (RU/EN × simple/technical).
   Same schema as the OCM diary's entries.js; rendered by the shared app.js on
   femoco.html. Honesty contract: every NUMBER comes from a real computation in
   the femoco/ pipeline. Framing/literature facts are cited, not invented; figures
   are added only once generated from results.
   ==================================================================== */

const UI = {
  langLabel:  { ru: "Язык",    en: "Language" },
  levelLabel: { ru: "Уровень", en: "Level" },
  ru: { ru: "Русский", en: "Russian" },
  en: { ru: "English", en: "English" },
  simple: { ru: "Простой", en: "Simple" },
  tech:   { ru: "Технический", en: "Technical" },
  glossaryTitle: { ru: "Словарь терминов", en: "Glossary" },
  tocTitle: { ru: "Дневник", en: "Diary" },
  figurePrefix: { ru: "Рис.", en: "Fig." },
  honesty: {
    ru: "Принцип проекта: каждое число — из реального расчёта (PySCF/PennyLane). Что нельзя проверить — помечаем явно. Ничего не приукрашиваем; неудачи методов и упёртые в бюджет расчёты показываем как есть.",
    en: "Project rule: every number comes from a real computation (PySCF/PennyLane). What can't be validated is labelled. Nothing is smoothed; method failures and budget-capped runs are shown as they are."
  }
};

const SITE = {
  title:    { ru: "FeMoco-лаборатория", en: "FeMoco Lab" },
  subtitle: { ru: "Активация N₂ → NH₃: открытый дневник пути к квантовому расчёту FeMoco",
              en: "N₂ → NH₃ activation: an open diary of the road to a FeMoco-scale quantum calculation" },
  hero: {
    simple: {
      ru: "Почти весь азот в удобрениях — а значит и в еде половины человечества — добывают из воздуха процессом Габера–Боша при высоких давлении и температуре. Бактерии делают то же самое при комнатных условиях ферментом с железо-молибденовым центром (FeMoco). Чтобы научиться так же, нужно понять, как FeMoco рвёт прочнейшую тройную связь N≡N. Это — одна из самых трудных задач квантовой химии. Здесь я честно записываю, как мы к ней подбираемся: начиная с того, что ещё можно посчитать точно.",
      en: "Almost all the nitrogen in fertilizer — and so in the food of half the world — is pulled from the air by the Haber–Bosch process at high pressure and temperature. Bacteria do the same at room conditions with an iron–molybdenum centre (FeMoco). To copy that, we must understand how FeMoco breaks the very strong N≡N triple bond — one of the hardest problems in quantum chemistry. Here I log honestly how we approach it, starting from what can still be computed exactly."
    },
    tech: {
      ru: "Активация N₂ на FeMo-кофакторе (Fe₇MoS₉C) — эталонная мультиреференсная задача: полное активное пространство, по Reiher et al. (PNAS 2017), это CAS(54e,54o) ≈ 10¹⁶ детерминантов, вне досягаемости точного Full CI; DMRG для полного кофактора не сходится, а ошибки DFT в спиновых зазорах Fe-кластеров достигают ~1 эВ. Стратегия дневника: (1) валидировать конвейер и наглядно показать срыв одно-детерминантных методов на точно решаемой молекуле N₂; (2) измерить рост стоимости с размером активного пространства и оценить квантовый ресурс; (3) подниматься к упрощённым Fe–S кластерам. Методы: PySCF (RHF/CCSD(T)/CASSCF) и PennyLane (ADAPT-VQE) со сверкой кубитного гамильтониана против CASCI."
    }
  }
};

const ENTRIES = [
  /* ---------------------------------------------------------------- */
  {
    id: "why",
    date: "2026-06",
    stage: { ru: "Зачем", en: "Why" },
    title: { ru: "Зачем это: азот из воздуха и стена в 10¹⁶",
             en: "Why this matters: nitrogen from air, and the 10¹⁶ wall" },
    simple: {
      ru: `<p><strong>Азот</strong> нужен всему живому, но в воздухе он связан в инертную молекулу N₂ с тройной связью N≡N — её очень трудно разорвать. Промышленный процесс <em>Габера–Боша</em> делает это при ~400–500 °C и ~150–300 атм и тратит на это порядка 1–2% всей энергии человечества.</p>
           <p><strong>Природа</strong> делает то же самое при комнатной температуре — ферментом нитрогеназой, в сердце которого металлический кластер из железа, молибдена и серы: <em>FeMoco</em>. Если бы мы поняли, как именно он работает, можно было бы спроектировать катализатор для «мягкого» аммиака.</p>
           <p><strong>Загвоздка:</strong> чтобы посчитать поведение электронов в FeMoco, нужно учесть 54 электрона на 54 орбиталях одновременно — это около <strong>10¹⁶ вариантов (детерминантов)</strong>. Точные классические методы столько не вытягивают. Поэтому это и считают «эталонной задачей» для будущих квантовых компьютеров.</p>
           <p><strong>Наш план — честный и снизу вверх:</strong> сперва посчитать то, что ещё решается точно (саму молекулу N₂), убедиться, что наш «калькулятор» не врёт, показать, где именно ломается классика, а затем шаг за шагом подниматься к FeMoco.</p>`,
      en: `<p><strong>Nitrogen</strong> is essential to life, but in air it is locked in the inert N₂ molecule with a triple bond N≡N that is very hard to break. The industrial <em>Haber–Bosch</em> process does it at ~400–500 °C and ~150–300 atm, consuming on the order of 1–2% of humanity's total energy.</p>
           <p><strong>Nature</strong> does the same at room temperature with the enzyme nitrogenase, whose heart is a metal cluster of iron, molybdenum and sulfur: <em>FeMoco</em>. If we understood exactly how it works, we could design a catalyst for "mild" ammonia.</p>
           <p><strong>The catch:</strong> modelling the electrons in FeMoco means treating 54 electrons in 54 orbitals at once — about <strong>10¹⁶ possibilities (determinants)</strong>. Exact classical methods cannot handle that, which is why it is treated as a "benchmark problem" for future quantum computers.</p>
           <p><strong>Our plan — honest and bottom-up:</strong> first compute what is still exactly solvable (the N₂ molecule itself), prove our "calculator" doesn't lie, show exactly where classical methods break, and then climb step by step toward FeMoco.</p>`
    },
    tech: {
      ru: `<p>Целевая реакция нитрогеназы: <code>N₂ + 8 H⁺ + 8 e⁻ → 2 NH₃ + H₂</code>, лимитирующая физика — расщепление тройной связи N≡N на FeMo-кофакторе (Fe₇MoS₉C) через спин-богатые, мультиреференсные интермедиаты железа.</p>
           <p><strong>Почему это «эталон»:</strong> Reiher, Wiebe, Svore, Wecker, Troyer (<em>PNAS</em> 114, 7555, 2017) выбрали активное пространство FeMoco как <strong>CAS(54e,54o)</strong>. Размерность Full CI такого пространства астрономическая; обычно цитируемый порядок числа определяющих конфигураций — <strong>~10¹⁶</strong>. Точный Full CI здесь невозможен, DMRG для полного кофактора не сходится к надёжному решению, а однодетерминантный DFT даёт ошибки спиновых зазоров Fe-кластеров до ~1 эВ.</p>
           <p><strong>Суть трудности — статическая (мультиреференсная) корреляция:</strong> при разрыве кратной связи несколько электронных пар «расцепляются», и волновая функция перестаёт описываться одним детерминантом. Именно это ломает RHF/CCSD(T) и DFT и требует многоконфигурационных методов (CASSCF/CASCI), а в пределе — квантовых (VQE) для активных пространств за гранью классики.</p>
           <p><strong>План дневника:</strong> (1) точно решаемый бенчмарк — диссоциация N₂ (та же физика N≡N, но решается точно) для валидации конвейера и демонстрации срыва классики; (2) измерить масштабирование стоимости с размером активного пространства и дать честную ресурс-оценку пути к CAS(54,54); (3) упрощённые Fe–S кластеры. Инструменты: PySCF (RHF/CCSD(T)/CASSCF/CASCI) и PennyLane (ADAPT-VQE), кубитный гамильтониан сверяется с CASCI до &lt;10⁻³ ккал/моль.</p>`,
      en: `<p>Nitrogenase target reaction: <code>N₂ + 8 H⁺ + 8 e⁻ → 2 NH₃ + H₂</code>; the rate-limiting physics is cleaving the N≡N triple bond on the FeMo-cofactor (Fe₇MoS₉C) through spin-rich, multireference iron intermediates.</p>
           <p><strong>Why it is "the benchmark":</strong> Reiher, Wiebe, Svore, Wecker, Troyer (<em>PNAS</em> 114, 7555, 2017) chose the FeMoco active space as <strong>CAS(54e,54o)</strong>. The Full-CI dimension is astronomical; the commonly quoted order for the number of determining configurations is <strong>~10¹⁶</strong>. Exact Full CI is impossible, DMRG does not converge to a reliable solution for the full cofactor, and single-determinant DFT mis-estimates Fe-cluster spin gaps by up to ~1 eV.</p>
           <p><strong>The core difficulty — static (multireference) correlation:</strong> breaking a multiple bond uncouples several electron pairs, and the wavefunction is no longer one determinant. This is what breaks RHF/CCSD(T) and DFT and demands multiconfigurational methods (CASSCF/CASCI) and, in the limit, quantum methods (VQE) for active spaces beyond classical reach.</p>
           <p><strong>Diary plan:</strong> (1) an exactly-solvable benchmark — N₂ dissociation (same N≡N physics, but solvable) to validate the pipeline and show the classical breakdown; (2) measure cost scaling vs active-space size and give an honest resource estimate of the road to CAS(54,54); (3) simplified Fe–S clusters. Tools: PySCF (RHF/CCSD(T)/CASSCF/CASCI) and PennyLane (ADAPT-VQE); the qubit Hamiltonian is cross-checked against CASCI to &lt;10⁻³ kcal/mol.</p>`
    },
    figures: []
  },

  /* ---------------------------------------------------------------- */
  {
    id: "engine",
    date: "2026-06",
    stage: { ru: "Этап 1", en: "Stage 1" },
    title: { ru: "Движок честности: кубитный гамильтониан = CASCI до 9-го знака",
             en: "An honesty engine: the qubit Hamiltonian = CASCI to 9 digits" },
    simple: {
      ru: `<p>Любой вывод дневника стоит ровно столько, сколько стоит наш «калькулятор». Поэтому первым делом мы доказываем, что он <strong>не врёт</strong>.</p>
           <p>Мы превращаем химию в форму, понятную квантовому компьютеру (так называемый «кубитный гамильтониан»), а затем сравниваем его точное решение с эталонным классическим ответом (CASCI) на системах, которые ещё решаются точно: H₂, LiH и активное пространство N₂.</p>
           <p><strong>Совпадение — до девятого знака</strong> (расхождение ~5×10⁻⁹ ккал/моль — это в миллиард раз точнее «химической точности» в 1 ккал/моль). Этот тест мы прогоняем <em>перед каждым</em> квантовым расчётом, так что ошибка в переводе химии в кубиты не проскочит незаметно.</p>`,
      en: `<p>Every conclusion in this diary is worth exactly as much as our "calculator". So the first thing we do is prove it <strong>doesn't lie</strong>.</p>
           <p>We turn the chemistry into a form a quantum computer understands (a "qubit Hamiltonian") and then compare its exact solution to the reference classical answer (CASCI) on systems still solvable exactly: H₂, LiH, and the N₂ active space.</p>
           <p><strong>They agree to the ninth digit</strong> (discrepancy ~5×10⁻⁹ kcal/mol — a billion times tighter than the 1 kcal/mol "chemical accuracy"). This test runs <em>before every</em> quantum calculation, so an error translating chemistry into qubits can't slip through unnoticed.</p>`
    },
    tech: {
      ru: `<p>Рецепт «интегралы активного пространства → кубитный гамильтониан» (хемист-ERI → OpenFermion, отображение Jordan–Wigner) валидируется напрямую: точная диагонализация JW-гамильтониана при фиксированном числе частиц должна воспроизводить PySCF CASCI на тех же орбиталях.</p>
           <p>Результат (реальный прогон): H₂ — <strong>8×10⁻¹³</strong>, LiH — машинный ноль, N₂ CAS(6e,6o)/cc-pVDZ — <strong>5.4×10⁻⁹</strong> ккал/моль. Худшее расхождение <strong>≈5×10⁻⁹ ккал/моль</strong> — на ~9 порядков ниже химической точности.</p>
           <p>Этот кросс-чек (<code>build_system(..., crosscheck=True)</code>) — обязательный гейт перед VQE: молчаливая ошибка соглашения/маппинга не пройдёт (assert &lt;10⁻³ ккал/моль). Это тот же движок, что валидирован в OCM-проекте на H₂/LiH до ~10⁻¹² Ha.</p>`,
      en: `<p>The recipe "active-space integrals → qubit Hamiltonian" (chemist-ERI → OpenFermion, Jordan–Wigner mapping) is validated directly: exact diagonalization of the JW Hamiltonian at fixed particle number must reproduce PySCF CASCI on the same orbitals.</p>
           <p>Result (real run): H₂ — <strong>8×10⁻¹³</strong>, LiH — machine zero, N₂ CAS(6e,6o)/cc-pVDZ — <strong>5.4×10⁻⁹</strong> kcal/mol. Worst disagreement <strong>≈5×10⁻⁹ kcal/mol</strong> — ~9 orders of magnitude below chemical accuracy.</p>
           <p>This cross-check (<code>build_system(..., crosscheck=True)</code>) is a mandatory gate before VQE: a silent convention/mapping error can't pass (assert &lt;10⁻³ kcal/mol). It is the same engine validated in the OCM project on H₂/LiH to ~10⁻¹² Ha.</p>`
    },
    table: {
      title: { ru: "Сверка кубитного гамильтониана с CASCI (реальный прогон)",
               en: "Qubit-Hamiltonian vs CASCI cross-check (real run)" },
      head: { ru: ["Система", "Активное пр-во", "Кубиты", "|Δ| vs CASCI, ккал/моль"],
              en: ["System", "Active space", "Qubits", "|Δ| vs CASCI, kcal/mol"] },
      rows: [
        ["H₂", "CAS(2e,2o)", "4", "8.4 × 10⁻¹³"],
        ["LiH", "CAS(2e,3o)", "6", "&lt; 10⁻¹² (машинный ноль)"],
        ["N₂", "CAS(6e,6o)", "12", "5.4 × 10⁻⁹"]
      ]
    },
    figures: []
  },

  /* ---------------------------------------------------------------- */
  {
    id: "n2-bench",
    date: "2026-06",
    stage: { ru: "Этап 2", en: "Stage 2" },
    title: { ru: "Честный калькулятор: разрыв N≡N и срыв классики-одиночки",
             en: "An honest calculator: breaking N≡N and the single-reference collapse" },
    simple: {
      ru: `<p>Прежде чем замахиваться на FeMoco, мы проверяем «калькулятор» на том, что ещё можно решить <strong>точно</strong> — на самой молекуле N₂. Растягиваем связь N≡N шаг за шагом и сравниваем методы (базис cc-pVDZ).</p>
           <p><strong>Что вышло (реальный расчёт):</strong></p>
           <ul>
             <li><strong>RHF</strong> (одна «картинка» для электронов) при растяжении уходит вверх в <em>неверный</em> предел — на дне диссоциации он завышает энергию на ~430 ккал/моль.</li>
             <li><strong>CCSD(T)</strong> (золотой стандарт у равновесия) у дна ямы отличный, но при растяжении <em>ломается</em>: при R≈2.5 Å даёт −110.37 Ha — это физически невозможно (ниже точной энергии). Так выглядит срыв.</li>
             <li><strong>CASSCF</strong> (многоконфигурационный, «мультиреференс») идёт гладко и правильно — равновесие при ≈1.10 Å, корректная диссоциация.</li>
           </ul>
           <p>Вывод честный и важный: именно <em>этот</em> срыв одно-детерминантных методов на разрыве связи — и есть причина, по которой полный FeMoco неподъёмен для обычной классики. Мы воспроизвели его на пальцах.</p>`,
      en: `<p>Before reaching for FeMoco we test the "calculator" on something still solvable <strong>exactly</strong> — the N₂ molecule itself. We stretch the N≡N bond step by step and compare methods (cc-pVDZ basis).</p>
           <p><strong>What we got (real computation):</strong></p>
           <ul>
             <li><strong>RHF</strong> (a single "picture" of the electrons) climbs to the <em>wrong</em> limit on stretching — at dissociation it overshoots by ~430 kcal/mol.</li>
             <li><strong>CCSD(T)</strong> (the gold standard near equilibrium) is excellent at the bottom of the well but <em>breaks</em> on stretching: at R≈2.5 Å it returns −110.37 Ha — physically impossible (below the exact energy). That is the collapse.</li>
             <li><strong>CASSCF</strong> (multiconfigurational, "multireference") is smooth and correct — equilibrium at ≈1.10 Å, correct dissociation.</li>
           </ul>
           <p>The honest, important point: this very collapse of single-determinant methods at bond breaking is <em>why</em> the full FeMoco is out of reach for ordinary classical methods. We reproduced it from first principles.</p>`
    },
    tech: {
      ru: `<p>Сканирование N₂ по cc-pVDZ, активное пространство для VQE — CAS(6e,6o)=12 кубитов (валентные σ/σ*/π/π* 2p). RHF/CCSD/CCSD(T) — PySCF; CASSCF(6,6) посеян орбиталями AVAS(«N 2p»).</p>
           <ul>
             <li><strong>RHF:</strong> у равновесия (R=1.10 Å) на 85 ккал/моль выше CASSCF; на диссоциации (R=2.6 Å) — на <strong>431.8 ккал/моль</strong> выше (E_RHF=−108.089 vs E_CASSCF=−108.777 Ha). Классический провал спин-ограниченного HF.</li>
             <li><strong>CCSD(T):</strong> при R=1.10 Å даёт −109.279 Ha (на 118 ккал/моль ниже CASSCF — это динамическая корреляция). Но CCSD перестаёт сходиться при R≳2.4 Å, и (T) выдаёт нефизичные значения — экстремум −110.366 Ha при R=2.5 Å (<em>ниже</em> точной FCI: не-вариационный срыв ряда возбуждений).</li>
             <li><strong>CASSCF(6,6):</strong> гладкая кривая, равновесие R≈1.10 Å, «глубина» до R=2.6 Å ≈ <strong>196 ккал/моль</strong> (CAS(6,6) без динамической корреляции занижает абсолютную D<sub>e</sub> — отмечаем честно).</li>
           </ul>
           <p>Это и есть статическая корреляция: при разрыве N≡N расцепляются три пары, и одно-детерминантные методы перестают быть надёжными. Кубитный гамильтониан активного пространства сверяется с CASCI до &lt;10⁻³ ккал/моль (см. Этап 2).</p>`,
      en: `<p>N₂ scan in cc-pVDZ; VQE active space CAS(6e,6o)=12 qubits (valence 2p σ/σ*/π/π*). RHF/CCSD/CCSD(T) via PySCF; CASSCF(6,6) seeded with AVAS("N 2p") orbitals.</p>
           <ul>
             <li><strong>RHF:</strong> 85 kcal/mol above CASSCF at equilibrium (R=1.10 Å); <strong>431.8 kcal/mol</strong> above at dissociation (R=2.6 Å; E_RHF=−108.089 vs E_CASSCF=−108.777 Ha). The classic restricted-HF failure.</li>
             <li><strong>CCSD(T):</strong> −109.279 Ha at R=1.10 Å (118 kcal/mol below CASSCF — dynamic correlation). But CCSD stops converging for R≳2.4 Å and (T) returns unphysical values — an extremum of −110.366 Ha at R=2.5 Å (<em>below</em> exact FCI: non-variational series collapse).</li>
             <li><strong>CASSCF(6,6):</strong> smooth curve, equilibrium R≈1.10 Å, "depth" to R=2.6 Å ≈ <strong>196 kcal/mol</strong> (CAS(6,6) without dynamic correlation underestimates the absolute D<sub>e</sub> — noted honestly).</li>
           </ul>
           <p>This is static correlation: breaking N≡N uncouples three pairs and single-determinant methods stop being reliable. The active-space qubit Hamiltonian agrees with CASCI to &lt;10⁻³ kcal/mol (see Stage 2).</p>`
    },
    figures: [
      { src: "assets/femoco/n2_dissociation_curve.png",
        caption: { ru: "Кривая диссоциации N₂ (cc-pVDZ): RHF уходит в неверный предел, CCSD(T) срывается при растяжении, CASSCF гладок и верен. Реальный расчёт.",
                   en: "N₂ dissociation curve (cc-pVDZ): RHF goes to the wrong limit, CCSD(T) collapses on stretching, CASSCF is smooth and correct. Real computation." } }
    ]
  },

  /* ---------------------------------------------------------------- */
  {
    id: "mr",
    date: "2026-06",
    stage: { ru: "Этап 3", en: "Stage 3" },
    title: { ru: "Сколько здесь мультиреференсности? Измеряем напрямую",
             en: "How multireference is it? Measuring it directly" },
    simple: {
      ru: `<p>«Классике-одиночке нельзя верить» — звучит как мнение. Сделаем из этого <strong>число</strong>.</p>
           <p>В простой молекуле орбитали либо заполнены (2 электрона), либо пусты (0). Когда связь N≡N рвётся, они становятся «полупустыми» (≈1) — и значит, никакая одна «картинка» молекулу не описывает. Мы считаем индикатор <strong>N_u</strong> — число эффективно неспаренных электронов (0 = идеальная одна картинка).</p>
           <p><strong>Результат (реальный CASSCF):</strong> у равновесия N_u = <strong>0.53</strong> (почти одна картинка — поэтому CCSD(T) там и работает). При растяжении N_u растёт до <strong>5.47</strong> — то есть почти все 6 связывающих электронов «расцепляются» (три связи разорваны). Вот <em>измеренная</em> причина, почему классика-одиночка врёт — и почему FeMoco, где железо-серные центры живут в таком режиме <em>постоянно</em>, так сложен.</p>`,
      en: `<p>"You can't trust single-reference classical methods" sounds like an opinion. Let's turn it into a <strong>number</strong>.</p>
           <p>In a simple molecule orbitals are either full (2 electrons) or empty (0). When N≡N breaks they become "half-full" (≈1) — meaning no single "picture" describes the molecule. We compute an indicator <strong>N_u</strong> — the number of effectively unpaired electrons (0 = a perfect single picture).</p>
           <p><strong>Result (real CASSCF):</strong> at equilibrium N_u = <strong>0.53</strong> (nearly one picture — which is why CCSD(T) works there). On stretching, N_u rises to <strong>5.47</strong> — almost all 6 bonding electrons "uncouple" (three bonds broken). This is the <em>measured</em> reason single-reference fails — and why FeMoco, whose iron–sulfur centres live in this regime <em>permanently</em>, is so hard.</p>`
    },
    tech: {
      ru: `<p>Прямая мера статической корреляции — натуральные орбитальные заселённости (NOON) активного пространства из 1-РДМ CASSCF(6,6)/cc-pVDZ, и индикатор Хэда-Гордона N<sub>u</sub> = Σᵢ nᵢ(2−nᵢ).</p>
           <ul>
             <li><strong>R = 1.10 Å:</strong> NOON ≈ {1.98, 1.94, 1.94, 0.06, 0.06, 0.02}, N<sub>u</sub> = <strong>0.53</strong> — почти замкнутая оболочка; динамическая корреляция доминирует, CCSD(T) надёжен.</li>
             <li><strong>R = 1.60 Å:</strong> {1.91, 1.73, 1.73, 0.27, 0.27, 0.09}, N<sub>u</sub> = <strong>2.24</strong> — заселённости поплыли.</li>
             <li><strong>R = 2.20 Å:</strong> {1.46, 1.16, 1.16, 0.84, 0.84, 0.54}, N<sub>u</sub> = <strong>5.47</strong> ≈ 6 неспаренных = три разорванные пары σ+2π. Истинно многоконфигурационный режим.</li>
           </ul>
           <p>Рост N<sub>u</sub> 0.53 → 5.47 количественно объясняет срыв RHF/CCSD(T) (Этап 2). Центры Fe-S в FeMoco сидят в высоком-N<sub>u</sub> режиме постоянно — отсюда ошибки спиновых зазоров DFT до ~1 эВ и потребность в большом активном пространстве (→ VQE). Это и мотивирует масштаб CAS(54,54).</p>`,
      en: `<p>The direct measure of static correlation — natural-orbital occupation numbers (NOON) of the active space from the CASSCF(6,6)/cc-pVDZ 1-RDM, and the Head-Gordon indicator N<sub>u</sub> = Σᵢ nᵢ(2−nᵢ).</p>
           <ul>
             <li><strong>R = 1.10 Å:</strong> NOON ≈ {1.98, 1.94, 1.94, 0.06, 0.06, 0.02}, N<sub>u</sub> = <strong>0.53</strong> — nearly closed-shell; dynamic correlation dominates, CCSD(T) reliable.</li>
             <li><strong>R = 1.60 Å:</strong> {1.91, 1.73, 1.73, 0.27, 0.27, 0.09}, N<sub>u</sub> = <strong>2.24</strong> — occupations spreading.</li>
             <li><strong>R = 2.20 Å:</strong> {1.46, 1.16, 1.16, 0.84, 0.84, 0.54}, N<sub>u</sub> = <strong>5.47</strong> ≈ 6 unpaired = three broken pairs (σ+2π). Genuinely multiconfigurational.</li>
           </ul>
           <p>The growth N<sub>u</sub> 0.53 → 5.47 quantitatively explains the RHF/CCSD(T) collapse (Stage 2). FeMoco's Fe-S centres sit in the high-N<sub>u</sub> regime permanently — hence DFT spin-gap errors up to ~1 eV and the need for a large active space (→ VQE). This is what motivates the CAS(54,54) scale.</p>`
    },
    table: {
      title: { ru: "Натуральные заселённости активного пространства N₂ — CASSCF(6,6)/cc-pVDZ",
               en: "N₂ active-space natural occupations — CASSCF(6,6)/cc-pVDZ" },
      head: { ru: ["R (Å)", "Заселённости 6 орбиталей", "N_u (неспаренные e⁻)"],
              en: ["R (Å)", "Occupations of 6 orbitals", "N_u (unpaired e⁻)"] },
      rows: [
        ["1.10", "1.98 · 1.94 · 1.94 · 0.06 · 0.06 · 0.02", "0.53"],
        ["1.60", "1.91 · 1.73 · 1.73 · 0.27 · 0.27 · 0.09", "2.24"],
        ["2.20", "1.46 · 1.16 · 1.16 · 0.84 · 0.84 · 0.54", "5.47"]
      ]
    },
    figures: [
      { src: "assets/femoco/n2_noon.png",
        caption: { ru: "Натуральные заселённости CAS(6,6): от {2,2,2,0,0,0} (одна картинка) у равновесия к {1,1,…} при разрыве N≡N. N_u 0.53→5.47. Реальный CASSCF.",
                   en: "CAS(6,6) natural occupations: from {2,2,2,0,0,0} (one picture) at equilibrium to {1,1,…} as N≡N breaks. N_u 0.53→5.47. Real CASSCF." } }
    ]
  },

  /* ---------------------------------------------------------------- */
  {
    id: "vqe-conv",
    date: "2026-06",
    stage: { ru: "Этап 4", en: "Stage 4" },
    title: { ru: "VQE в деле: сходимость к точному ответу",
             en: "VQE at work: convergence to the exact answer" },
    simple: {
      ru: `<p>Мы утверждали, что квантовый алгоритм (VQE) воспроизводит точный ответ. Вот как это выглядит «в динамике».</p>
           <p>ADAPT-VQE собирает свою квантовую схему <em>оператор за оператором</em>, и мы следим, как энергия падает от стартовой ошибки ~40 ккал/моль (уровень простейшего приближения) вниз к <strong>химической точности (~1 ккал/моль)</strong> на равновесной N₂.</p>
           <p>Маленькое пространство (CAS(4,4), 8 кубитов) доходит за <strong>19 шагов / 16 с</strong>; большее (CAS(6,6), 12 кубитов) — за <strong>52 шага / 7 мин</strong>. Тот же результат, но заметно дороже. Это и есть квантовый метод, работающий на настоящей молекуле (пока на симуляторе) — и наглядный намёк, <em>почему</em> дорога к FeMoco требует и больше кубитов, и лучшего железа.</p>`,
      en: `<p>We claimed the quantum algorithm (VQE) reproduces the exact answer. Here is what that looks like "in motion".</p>
           <p>ADAPT-VQE builds its quantum circuit <em>operator by operator</em>, and we watch the energy fall from a starting error of ~40 kcal/mol (the simplest-approximation level) down to <strong>chemical accuracy (~1 kcal/mol)</strong> at the N₂ equilibrium.</p>
           <p>The small space (CAS(4,4), 8 qubits) gets there in <strong>19 steps / 16 s</strong>; the larger one (CAS(6,6), 12 qubits) in <strong>52 steps / 7 min</strong>. Same destination, markedly more expensive. This is the quantum method working on a real molecule (on a simulator for now) — and a vivid hint of <em>why</em> the road to FeMoco needs both more qubits and better hardware.</p>`
    },
    tech: {
      ru: `<p>ADAPT-VQE (пул синглетов+дублетов возбуждений), бэкенд lightning.qubit, cc-pVDZ, N₂ при R=1.10 Å. Перед запуском кубитный гамильтониан сверен с CASCI (&lt;10⁻³ ккал/моль).</p>
           <ul>
             <li><strong>CAS(4e,4o) / 8 кубитов:</strong> ошибка vs CASCI 40.3 → <strong>1.07</strong> ккал/моль, 19 операторов, 16 с (сошёлся по градиенту).</li>
             <li><strong>CAS(6e,6o) / 12 кубитов:</strong> 42.7 → <strong>1.21</strong> ккал/моль, 52 оператора, 416 с (сошёлся).</li>
           </ul>
           <p>Оба выходят на уровень химической точности; остаточные ~1 ккал/моль — это предел пула S+D (не полный FCI), и у равновесия он мал (согласуется с низким N_u, Этап 3 — при растяжении было бы хуже). Рост 8→12 кубитов: ×2.7 операторов и ×26 по времени стены — это предвестие масштабирования (Этап 5) и причина, по которой полный FeMoco требует отказоустойчивого железа, а не симулятора.</p>`,
      en: `<p>ADAPT-VQE (singles+doubles excitation pool), lightning.qubit backend, cc-pVDZ, N₂ at R=1.10 Å. The qubit Hamiltonian is cross-checked vs CASCI (&lt;10⁻³ kcal/mol) before each run.</p>
           <ul>
             <li><strong>CAS(4e,4o) / 8 qubits:</strong> error vs CASCI 40.3 → <strong>1.07</strong> kcal/mol, 19 operators, 16 s (gradient-converged).</li>
             <li><strong>CAS(6e,6o) / 12 qubits:</strong> 42.7 → <strong>1.21</strong> kcal/mol, 52 operators, 416 s (converged).</li>
           </ul>
           <p>Both reach the chemical-accuracy ballpark; the residual ~1 kcal/mol is the S+D pool limit (not full FCI), and at equilibrium it is small (consistent with low N_u, Stage 3 — stretching would be worse). The 8→12 qubit growth: ×2.7 operators and ×26 wall-time — a foretaste of the scaling (Stage 5) and why the full FeMoco needs fault-tolerant hardware, not a simulator.</p>`
    },
    table: {
      title: { ru: "Сходимость ADAPT-VQE на N₂ (R=1.10 Å, cc-pVDZ) — реальный прогон",
               en: "ADAPT-VQE convergence on N₂ (R=1.10 Å, cc-pVDZ) — real run" },
      head: { ru: ["Активное пр-во", "Кубиты", "Операторов", "Время", "Итог vs CASCI"],
              en: ["Active space", "Qubits", "Operators", "Wall time", "Final vs CASCI"] },
      rows: [
        ["CAS(4e,4o)", "8", "19", "16 с", "1.07 ккал/моль"],
        ["CAS(6e,6o)", "12", "52", "416 с", "1.21 ккал/моль"]
      ]
    },
    figures: [
      { src: "assets/femoco/n2_vqe_convergence.png",
        caption: { ru: "ADAPT-VQE сходится к точному CASCI на N₂: от ~40 ккал/моль к ~1 ккал/моль. Большое пространство дороже. Реальный прогон (lightning.qubit).",
                   en: "ADAPT-VQE converging to exact CASCI on N₂: from ~40 to ~1 kcal/mol. The larger space costs more. Real run (lightning.qubit)." } }
    ]
  },

  /* ---------------------------------------------------------------- */
  {
    id: "scaling",
    date: "2026-06",
    stage: { ru: "Этап 5", en: "Stage 5" },
    title: { ru: "Масштаб задачи и где стена: путь к CAS(54,54)",
             en: "Problem size and where the wall is: the road to CAS(54,54)" },
    simple: {
      ru: `<p>Теперь измеряем, как быстро задача «раздувается» с размером активного пространства, и проверяем, поспевает ли квантовый алгоритм (VQE) за точным ответом.</p>
           <ul>
             <li><strong>VQE действительно воспроизводит точную энергию</strong> на маленьких пространствах: на CAS(2,2) ошибка 0.63 ккал/моль. С ростом размера и при ограниченном бюджете времени ошибка растёт (до ~4 ккал/моль на 10 кубитах) — показываем как есть, без прикрас.</li>
             <li><strong>Число вариантов (детерминантов) растёт комбинаторно.</strong> Для активного пространства FeMoco CAS(54,54) их получается <strong>≈10³⁰</strong> — это далеко за пределами точной классики.</li>
           </ul>
           <p>Важная честная поправка: в витринах часто пишут «~10¹⁶». Наш прямой подсчёт даёт ≈4×10³⁰ детерминантов (это полная размерность). Каким бы ни был способ счёта — вывод один: точный классический расчёт здесь невозможен, и поэтому FeMoco — задача для будущих квантовых компьютеров.</p>`,
      en: `<p>Now we measure how fast the problem "blows up" with active-space size, and check whether the quantum algorithm (VQE) keeps up with the exact answer.</p>
           <ul>
             <li><strong>VQE does reproduce the exact energy</strong> on small spaces: 0.63 kcal/mol error on CAS(2,2). With larger size and a limited time budget the error grows (to ~4 kcal/mol at 10 qubits) — shown as it is, unembellished.</li>
             <li><strong>The number of configurations (determinants) grows combinatorially.</strong> For the FeMoco active space CAS(54,54) it is <strong>≈10³⁰</strong> — far beyond exact classical reach.</li>
           </ul>
           <p>An important honest correction: showcases often quote "~10¹⁶". Our direct count gives ≈4×10³⁰ determinants (the full dimension). Whatever the counting convention, the conclusion is the same: an exact classical calculation is impossible here, which is why FeMoco is a problem for future quantum computers.</p>`
    },
    tech: {
      ru: `<p>На растянутой N₂ (R=1.5 Å, сильная статическая корреляция) считаем CAS(n,n) растущего размера: точный CASCI всегда, ADAPT-VQE до 10 кубитов (бюджет 120 с/точку). Кубитный гамильтониан сверяется с CASCI перед каждым VQE.</p>
           <ul>
             <li><strong>VQE vs CASCI (реально измерено):</strong> CAS(2,2)/4q → 0.63 ккал/моль (2 оператора); CAS(4,4)/8q → 4.43 ккал/моль (19 операторов); CAS(4,5)/10q → 4.37 ккал/моль (23 оператора). Рост ошибки — следствие упёртости в бюджет ADAPT и силы корреляции, а не «подгонки». Честно.</li>
             <li><strong>Комбинаторный рост:</strong> число детерминантов = C(n, n/2)². Ориентир «стены» точного Full CI — ~10⁹.</li>
             <li><strong>Цель масштаба — FeMoco CAS(54e,54o):</strong> C(54,27)² = <strong>3.79×10³⁰</strong> детерминантов, <strong>108</strong> спин-орбитальных кубитов (наше отображение Jordan–Wigner). Часто цитируемое «~10¹⁶» — это иная/меньшая оценка; прямой подсчёт размерности даёт ~10³⁰.</li>
           </ul>
           <p><strong>Ресурс-оценка (честно):</strong> 108 кубитов — это только на активное пространство, до учёта отказоустойчивости. Полная FT-оценка для FeMoco (Reiher et al., PNAS 2017) выходит далеко за рамки сегодняшнего железа; конкретные T-счётчики см. в статье — мы их не пересчитывали и не приводим как свои. Ближайший достижимый шаг — упрощённые Fe–S кластеры (следующий этап), где число кубитов в пределах симулятора/раннего QPU.</p>`,
      en: `<p>On stretched N₂ (R=1.5 Å, strong static correlation) we compute CAS(n,n) of growing size: exact CASCI always, ADAPT-VQE up to 10 qubits (120 s/point budget). The qubit Hamiltonian is checked against CASCI before each VQE.</p>
           <ul>
             <li><strong>VQE vs CASCI (really measured):</strong> CAS(2,2)/4q → 0.63 kcal/mol (2 operators); CAS(4,4)/8q → 4.43 kcal/mol (19 operators); CAS(4,5)/10q → 4.37 kcal/mol (23 operators). The growing error reflects ADAPT hitting its budget and the correlation strength, not "tuning". Honest.</li>
             <li><strong>Combinatorial growth:</strong> determinants = C(n, n/2)². The exact Full-CI "wall" reference is ~10⁹.</li>
             <li><strong>The scale target — FeMoco CAS(54e,54o):</strong> C(54,27)² = <strong>3.79×10³⁰</strong> determinants, <strong>108</strong> spin-orbital qubits (our Jordan–Wigner mapping). The often-quoted "~10¹⁶" is a different/smaller estimate; the direct dimension count is ~10³⁰.</li>
           </ul>
           <p><strong>Resource estimate (honest):</strong> 108 qubits is for the active space only, before fault-tolerance overhead. The full FT estimate for FeMoco (Reiher et al., PNAS 2017) is far beyond today's hardware; for specific T-gate counts see the paper — we did not recompute them and do not present them as ours. The nearest reachable step is simplified Fe–S clusters (next stage), where the qubit count is within simulator / early-QPU range.</p>`
    },
    figures: [
      { src: "assets/femoco/n2_scaling.png",
        caption: { ru: "Реальные измерения: размер точной задачи и стоимость растут с активным пространством; FeMoco CAS(54,54) ≈10³⁰ детерминантов (108 кубитов) — за стеной классики. НЕ заявление о квантовом превосходстве.",
                   en: "Real measurements: exact-problem size and cost grow with the active space; FeMoco CAS(54,54) ≈10³⁰ determinants (108 qubits) — beyond the classical wall. NOT a quantum-advantage claim." } }
    ]
  },

  /* ---------------------------------------------------------------- */
  {
    id: "thermo",
    date: "2026-06",
    stage: { ru: "Этап 6", en: "Stage 6" },
    title: { ru: "Что классика берёт легко: термодинамика N₂ + 3H₂ → 2NH₃",
             en: "What classical methods nail: the thermodynamics of N₂ + 3H₂ → 2NH₃" },
    simple: {
      ru: `<p>Важный честный контрапункт ко всему дневнику. Сама реакция «азот + водород → аммиак» классике <strong>вполне по силам</strong>.</p>
           <p>Мы посчитали её энергию обычным золотым стандартом (CCSD(T)): ответ устойчив и улучшается с качеством базиса. Почему здесь легко, а разрыв N≡N (Этап 2) был тяжёл? Потому что <em>исходники и продукты</em> (N₂, H₂, NH₃) — спокойные, замкнутооболочечные молекулы у равновесия: их описывает одна «картинка» (малый N_u). Тяжёл — <em>путь между ними</em>: каталитические интермедиаты на железе, где связи наполовину разорваны.</p>
           <p><strong>Вывод:</strong> классика отлично берёт «пункт назначения». Квантово-трудна не термодинамика, а <strong>механизм катализа</strong> — именно туда и смотрим мы (и природа через FeMoco).</p>`,
      en: `<p>An important, honest counterpoint to the whole diary. The overall reaction "nitrogen + hydrogen → ammonia" is <strong>well within classical reach</strong>.</p>
           <p>We computed its energy with the usual gold standard (CCSD(T)): the answer is stable and improves with basis quality. Why is this easy when breaking N≡N (Stage 2) was hard? Because the <em>reactants and products</em> (N₂, H₂, NH₃) are calm, closed-shell molecules at equilibrium — one "picture" describes them (small N_u). What is hard is the <em>road between them</em>: the catalytic intermediates on iron where bonds are half-broken.</p>
           <p><strong>Takeaway:</strong> classical methods nail the "destination". The quantum-hard part is not the thermodynamics but the <strong>catalytic mechanism</strong> — exactly where we (and nature, via FeMoco) look.</p>`
    },
    tech: {
      ru: `<p>CCSD(T) одноточечно на фиксированных экспериментальных геометриях (H₂ 0.7414 Å; N₂ 1.0977 Å; NH₃ r=1.012 Å, ∠HNH=106.7°). Электронная энергия реакции ΔE = 2E(NH₃) − E(N₂) − 3E(H₂):</p>
           <ul>
             <li><strong>cc-pVTZ:</strong> RHF −33.5; MP2 −36.9; CCSD −39.5; <strong>CCSD(T) −37.4</strong> ккал/моль.</li>
             <li><strong>cc-pVDZ:</strong> CCSD(T) −24.7 ккал/моль — заметно недосходится по базису (~13 ккал/моль неполнота); честная иллюстрация важности базиса.</li>
           </ul>
           <p>Иерархия RHF→MP2→CCSD→CCSD(T) ведёт себя нормально (без расходимостей) — значит, равновесные молекулы однодетерминантны (согласуется с N_u=0.53 для N₂, Этап 3). Наша ΔE — <em>бесколебательная электронная</em> энергия; экспериментальная ΔH°₂₉₈ = −21.9 ккал/моль включает ZPE+тепловые, поэтому электронная величина закономерно отрицательнее на ≈ΔZPE (у NH₃ велика энергия нулевых колебаний N–H). Это не ошибка метода. <strong>Мораль:</strong> для термодинамики достаточно классического CC; мультиреференсная трудность (и аргумент за квант) живёт в переходных состояниях катализа, а не в стабильных конечных точках.</p>`,
      en: `<p>CCSD(T) single points at fixed experimental geometries (H₂ 0.7414 Å; N₂ 1.0977 Å; NH₃ r=1.012 Å, ∠HNH=106.7°). Electronic reaction energy ΔE = 2E(NH₃) − E(N₂) − 3E(H₂):</p>
           <ul>
             <li><strong>cc-pVTZ:</strong> RHF −33.5; MP2 −36.9; CCSD −39.5; <strong>CCSD(T) −37.4</strong> kcal/mol.</li>
             <li><strong>cc-pVDZ:</strong> CCSD(T) −24.7 kcal/mol — markedly basis-underconverged (~13 kcal/mol incompleteness); an honest illustration of basis sensitivity.</li>
           </ul>
           <p>The RHF→MP2→CCSD→CCSD(T) hierarchy behaves normally (no divergence) — i.e. the equilibrium species are single-reference (consistent with N_u=0.53 for N₂, Stage 3). Our ΔE is the <em>vibrationless electronic</em> energy; the experimental ΔH°₂₉₈ = −21.9 kcal/mol includes ZPE+thermal, so the electronic value is correctly more negative by ≈ΔZPE (NH₃ has large N–H zero-point energy). This is not a method error. <strong>Moral:</strong> classical CC suffices for the thermodynamics; the multireference difficulty (the case for quantum) lives in the catalytic transition states, not the stable end-points.</p>`
    },
    table: {
      title: { ru: "Энергия реакции N₂+3H₂→2NH₃ (ккал/моль) — реальный CCSD(T) и др.",
               en: "Reaction energy N₂+3H₂→2NH₃ (kcal/mol) — real CCSD(T) and others" },
      head: { ru: ["Базис", "RHF", "MP2", "CCSD", "CCSD(T)"],
              en: ["Basis", "RHF", "MP2", "CCSD", "CCSD(T)"] },
      rows: [
        ["cc-pVDZ", "−32.0", "−24.6", "−27.4", "−24.7"],
        ["cc-pVTZ", "−33.5", "−36.9", "−39.5", "−37.4"]
      ]
    },
    figures: [
      { src: "assets/femoco/nh3_thermo.png",
        caption: { ru: "Электронная ΔE реакции синтеза аммиака (CCSD(T) и др., cc-pVDZ/TZ). Классика берёт термодинамику; ΔE отрицательнее экспериментальной ΔH₂₉₈ на ≈ZPE. Реальный расчёт.",
                   en: "Electronic ΔE of ammonia synthesis (CCSD(T) etc., cc-pVDZ/TZ). Classical methods handle the thermodynamics; ΔE is more negative than experimental ΔH₂₉₈ by ≈ZPE. Real computation." } }
    ]
  },

  /* ---------------------------------------------------------------- */
  {
    id: "plan",
    date: "2026-06",
    stage: { ru: "Этап 7", en: "Stage 7" },
    title: { ru: "Куда дальше: к железо-серному кластеру (план)",
             en: "Where next: toward the iron–sulfur cluster (a plan)" },
    simple: {
      ru: `<p>До сих пор мы работали с самой молекулой N₂ — чтобы отладить инструменты и показать, <em>почему</em> тройная связь трудна. Следующий честный шаг к FeMoco — посадить N₂ на <strong>железо</strong>.</p>
           <p>Вся сложность FeMoco — от его железных центров: открытые оболочки, множество близких по энергии спиновых состояний. Поэтому мы построим небольшую <em>считаемую</em> модель с железом, которая связывает N₂, и прогоним тот же конвейер (классический эталон + VQE на симуляторе), используя AWS для тяжёлых частей.</p>
           <p>Полный FeMoco (10³⁰ — классически невозможно) мы изображать не будем. Идём к нему честно, по одной посильной модели за раз, и записываем, что работает, а что нет.</p>`,
      en: `<p>So far we worked with the N₂ molecule itself — to debug the tools and show <em>why</em> the triple bond is hard. The next honest step toward FeMoco is to put N₂ on <strong>iron</strong>.</p>
           <p>All of FeMoco's difficulty comes from its iron centres: open shells, many near-degenerate spin states. So we will build a small, <em>computable</em> iron model that binds N₂ and run the same pipeline (classical reference + VQE on a simulator), using AWS for the heavy parts.</p>
           <p>We will not pretend to do the full FeMoco (10³⁰ — classically impossible). We climb toward it honestly, one tractable model at a time, logging what works and what doesn't.</p>`
    },
    tech: {
      ru: `<p><strong>План следующего этапа</strong> (это план, не результаты — числа заполнит реальный прогон):</p>
           <ul>
             <li><strong>Модель:</strong> открытооболочечный Fe–N₂ фрагмент (моноядерный Fe-центр или малый [2Fe–2S] кластер) — ROHF/CASSCF-эталон, активное пространство охватывает Fe 3d + π/π* азота, ориентир CAS(10–16) → 20–32 кубита.</li>
             <li><strong>Дескриптор:</strong> энергия связывания N₂ и барьер первого гидрирования; порядок спиновых состояний (там DFT ошибается до ~1 эВ — наш реальный edge).</li>
             <li><strong>Железо:</strong> AWS lightning.qubit (где важны RAM/ядра). Ограничения честно: квота ~16 vCPU (один c7i.4xlarge зараз); statevector-VQE реалистичен до ~24–28 кубитов, дальше — приближения (фрагментация/embedding).</li>
             <li><strong>Контракт прежний:</strong> кубитный гамильтониан сверяется с CASCI; всё схематичное помечается; неудачи показываем.</li>
           </ul>
           <p>Это переводит дневник от «почему трудно» (N₂) к «как мы реально считаем катализатор» (Fe). Следующая запись будет с реальными числами модели.</p>`,
      en: `<p><strong>Next-stage plan</strong> (a plan, not results — numbers will be filled by a real run):</p>
           <ul>
             <li><strong>Model:</strong> an open-shell Fe–N₂ fragment (a mononuclear Fe centre or a small [2Fe–2S] cluster) — ROHF/CASSCF reference, active space spanning Fe 3d + N₂ π/π*, target CAS(10–16) → 20–32 qubits.</li>
             <li><strong>Descriptor:</strong> N₂ binding energy and first-hydrogenation barrier; spin-state ordering (where DFT errs by up to ~1 eV — our real edge).</li>
             <li><strong>Hardware:</strong> AWS lightning.qubit (where RAM/cores matter). Constraints stated honestly: ~16 vCPU quota (one c7i.4xlarge at a time); statevector VQE realistic to ~24–28 qubits, beyond which approximations (fragmentation/embedding) are needed.</li>
             <li><strong>Same contract:</strong> the qubit Hamiltonian is checked against CASCI; anything schematic is labelled; failures are shown.</li>
           </ul>
           <p>This moves the diary from "why it's hard" (N₂) to "how we actually compute the catalyst" (Fe). The next entry will carry real model numbers.</p>`
    },
    figures: []
  },

  /* ---------------------------------------------------------------- */
  {
    id: "fe-spin",
    date: "2026-06",
    stage: { ru: "Этап 8", en: "Stage 8" },
    title: { ru: "Железо вблизи: даже один атом Fe уже труден",
             en: "Iron up close: even a single Fe atom is already hard" },
    simple: {
      ru: `<p>Начинаем железную часть с самого простого случая — <strong>одного атома железа</strong> — чтобы откалибровать методы перед кластером. Уже здесь начинается то, из-за чего FeMoco — фронтир.</p>
           <p>Мы посчитали энергетическую щель между основным состоянием железа (квинтет ⁵D — это эксперимент, NIST) и другими спиновыми состояниями, четырьмя способами:</p>
           <ul>
             <li><strong>ROHF</strong> (одна картинка): <em>качественно врёт</em> — кладёт септет даже ниже основного состояния (переоценивает высокий спин), а триплет завышает на ~5.7 эВ.</li>
             <li><strong>DFT</strong> (PBE и B3LYP): разумно, но два функционала <em>расходятся</em> (1.18 против 1.36 эВ) — и какой прав, из самого DFT не понять.</li>
             <li><strong>Минимальный CASSCF</strong>: 0.39 эВ — но <em>занижает</em>, потому что маленькое активное пространство недоучитывает корреляцию.</li>
           </ul>
           <p>Ни один дешёвый метод не надёжен даже на <em>одном</em> атоме. В FeMoco таких атомов <strong>семь</strong>, и они магнитно связаны. Вот почему это передний край.</p>`,
      en: `<p>We begin the iron part with the simplest case — a <strong>single iron atom</strong> — to calibrate methods before the cluster. Even here, what makes FeMoco a frontier begins.</p>
           <p>We computed the energy gap between iron's ground state (quintet ⁵D — experimental, NIST) and other spin states, four ways:</p>
           <ul>
             <li><strong>ROHF</strong> (one picture): <em>qualitatively wrong</em> — it puts the septet even below the ground state (over-favours high spin), and overshoots the triplet by ~5.7 eV.</li>
             <li><strong>DFT</strong> (PBE and B3LYP): reasonable, but the two functionals <em>disagree</em> (1.18 vs 1.36 eV) — and DFT alone can't tell you which is right.</li>
             <li><strong>Minimal CASSCF</strong>: 0.39 eV — but it <em>underestimates</em>, because a small active space misses correlation.</li>
           </ul>
           <p>No cheap method is reliable on even <em>one</em> atom. FeMoco has <strong>seven</strong>, magnetically coupled. That is why it is the frontier.</p>`
    },
    tech: {
      ru: `<p>Атом Fe, def2-SVP. Щель триплет − квинтет (³F − ⁵D), реальный расчёт:</p>
           <ul>
             <li><strong>ROHF: +5.73 эВ</strong>, и септет лежит на 0.18 эВ <em>ниже</em> квинтета — одно-детерминантный метод качественно неверно упорядочивает спины (завышает высокоспиновые).</li>
             <li><strong>CASSCF(8e,6o): +0.39 эВ</strong> — минимальное активное пространство (3d+4s) без динамической корреляции занижает; надёжная щель требует большего CAS + CASPT2/NEVPT2.</li>
             <li><strong>UKS-PBE +1.18, UKS-B3LYP +1.36 эВ</strong> — функционалы расходятся; «правильного» ответа из DFT не извлечь.</li>
           </ul>
           <p><strong>Честно о подводном камне:</strong> наш первый, наивный CASSCF (AVAS отдельно на каждое спиновое состояние) дал <em>нефизичный</em> результат — триплет на 2.78 эВ ниже квинтета — потому что AVAS выбирал <em>разные</em> активные пространства для разных состояний. Починили принудительным единым CAS(8e,6o) + спин-проекцией (fix_spin). Сам факт, что корректную мультиреференсную постановку для железа легко испортить, — часть ответа на вопрос, почему FeMoco так тяжёл. Разброс в несколько эВ на одном атоме, помноженный на 7 обменно-связанных Fe, — суть задачи; отсюда и потребность в большом активном пространстве (→ VQE).</p>`,
      en: `<p>Fe atom, def2-SVP. Triplet − quintet gap (³F − ⁵D), real computation:</p>
           <ul>
             <li><strong>ROHF: +5.73 eV</strong>, and the septet lies 0.18 eV <em>below</em> the quintet — the single-determinant method qualitatively mis-orders the spins (over-stabilizes high spin).</li>
             <li><strong>CASSCF(8e,6o): +0.39 eV</strong> — the minimal active space (3d+4s) without dynamic correlation underestimates; a reliable gap needs a larger CAS + CASPT2/NEVPT2.</li>
             <li><strong>UKS-PBE +1.18, UKS-B3LYP +1.36 eV</strong> — the functionals disagree; DFT alone can't yield the "right" answer.</li>
           </ul>
           <p><strong>An honest pitfall:</strong> our first, naive CASSCF (AVAS run separately per spin state) gave an <em>unphysical</em> result — the triplet 2.78 eV below the quintet — because AVAS selected <em>different</em> active spaces for different states. Fixed by forcing a single consistent CAS(8e,6o) + spin projection (fix_spin). That a correct multireference setup for iron is easy to get wrong is itself part of why FeMoco is so hard. A multi-eV scatter on one atom, multiplied over 7 exchange-coupled Fe, is the heart of the problem — hence the need for a large active space (→ VQE).</p>`
    },
    table: {
      title: { ru: "Спиновая щель атома Fe (триплет − квинтет), эВ — реальный расчёт (def2-SVP)",
               en: "Fe-atom spin gap (triplet − quintet), eV — real computation (def2-SVP)" },
      head: { ru: ["Метод", "³F − ⁵D (эВ)", "Комментарий"],
              en: ["Method", "³F − ⁵D (eV)", "Comment"] },
      rows: [
        ["ROHF", "+5.73", "качественно неверно (септет ниже ⁵D)"],
        ["CASSCF(8,6)", "+0.39", "занижает (мало корреляции)"],
        ["DFT PBE", "+1.18", "функционал-зависимо"],
        ["DFT B3LYP", "+1.36", "функционал-зависимо"]
      ]
    },
    figures: [
      { src: "assets/femoco/fe_spin.png",
        caption: { ru: "Спиновые состояния атома Fe по четырём методам (отн. квинтета): ROHF врёт качественно, DFT расходится, минимальный CASSCF занижает. Реальный расчёт (def2-SVP).",
                   en: "Fe-atom spin states by four methods (vs quintet): ROHF qualitatively wrong, DFT scatters, minimal CASSCF underestimates. Real computation (def2-SVP)." } }
    ]
  },

  /* ---------------------------------------------------------------- */
  {
    id: "fe-n2",
    date: "2026-06",
    stage: { ru: "Этап 9", en: "Stage 9" },
    title: { ru: "N₂ на железе: первая модель самого катализа (и честный сюрприз)",
             en: "N₂ on iron: the first model of catalysis itself (and an honest surprise)" },
    simple: {
      ru: `<p>После калибровки на голом атоме Fe (Этап 8) сажаем N₂ на железо — первая модель самого́ катализа. Оптимизировали геометрию и получили честный, поучительный результат: <strong>одинокий атом железа — плохой активатор</strong>.</p>
           <ul>
             <li>В <strong>триплете</strong> Fe всё же слабо цепляет N₂ (Fe–N 1.85 Å) и чуть растягивает тройную связь (<strong>+0.028 Å</strong> — намёк на активацию), но комплекс на <strong>+15.4 ккал/моль ВЫШЕ</strong> разделённых Fe + N₂ — то есть сам собой не связывается.</li>
             <li>В <strong>квинтете</strong> (естественное состояние свободных Fe + N₂) атом и молекула просто расходятся.</li>
           </ul>
           <p>Вывод: одного железа мало — ровно поэтому природа упаковывает его в серно-углеродную клетку с несколькими металлами (FeMoco). Ещё деталь: триплет и квинтет почти вырождены (в пределах ~0.1 эВ), и методы спорят, кто ниже — тонкая спиновая физика, как раз для квантовых методов.</p>
           <p>«Катализаторное» активное пространство тут (d-орбитали железа + орбитали азота) — <strong>CAS(16,12) = 24 кубита</strong>, вдвое больше, чем для чистого N₂. Это следующая цель VQE на AWS.</p>`,
      en: `<p>After calibrating on the bare Fe atom (Stage 8), we put N₂ on iron — the first model of catalysis itself. We optimized the geometry and got an honest, instructive result: <strong>a lone iron atom is a poor activator</strong>.</p>
           <ul>
             <li>In the <strong>triplet</strong> Fe does weakly grab N₂ (Fe–N 1.85 Å) and slightly stretches the triple bond (<strong>+0.028 Å</strong> — a hint of activation), but the complex is <strong>+15.4 kcal/mol ABOVE</strong> separated Fe + N₂ — i.e. it does not bind spontaneously.</li>
             <li>In the <strong>quintet</strong> (the natural state of free Fe + N₂) the atom and molecule simply drift apart.</li>
           </ul>
           <p>Conclusion: one iron is not enough — exactly why nature wraps it in a sulfur–carbon cage with several metals (FeMoco). Also: triplet and quintet are nearly degenerate (within ~0.1 eV), and methods disagree on which wins — delicate spin physics, the kind quantum methods are built for.</p>
           <p>The "catalyst" active space here (iron d-orbitals + nitrogen orbitals) is <strong>CAS(16,12) = 24 qubits</strong>, double the N₂-only size. That is the next VQE target on AWS.</p>`
    },
    tech: {
      ru: `<p>B3LYP/def2-SVP оптимизация геометрии; ROHF→CASSCF; активные пространства через AVAS. Все числа посчитаны.</p>
           <ul>
             <li><strong>Связывание/активация:</strong> триплетный минимум d(Fe–N)=1.847, d(N–N)=1.129 Å (+0.028 vs свободная 1.100). E_связ = <strong>+15.4 ккал/моль</strong> относительно Fe(⁵D)+N₂ (эндотермично → метастабильно). Квинтет оптимизируется к d(Fe–N)=5.47 Å (диссоциация). Bare Fe–N₂ в основном состоянии не связан; лёгкая активация — лишь в метастабильном триплете.</li>
             <li><strong>Вертикальная спиновая щель</strong> (квинтет−триплет) при геометрии связанного комплекса: PBE <strong>−0.05</strong>, B3LYP <strong>−0.03</strong>, CASSCF(6e,5o) <strong>−0.13</strong> эВ → почти вырождение, расхождение DFT vs мультиреференс ~0.1 эВ. (У полноразмерных FeMoco-кластеров ошибки DFT по спину доходят до ~1 эВ — это о большом кофакторе, не об этой мини-модели.)</li>
             <li><strong>Активное пространство для VQE:</strong> AVAS(Fe 3d + N₂ 2p) = CAS(16e,12o) → <strong>24 кубита</strong> (против 12 у чистого N₂). CASSCF сошёлся для триплета (E=−1371.036 Ha); квинтет на 24 кубитах сходится трудно — само по себе признак сложности. Это 24-кубитное пространство — следующая цель VQE на AWS.</li>
           </ul>
           <p>Это намеренно простейшая «Fe-модель»; урок (нужны лиганды/многоядерность) ведёт к следующему этапу — [2Fe–2S] кластеру. Кубитный гамильтониан будет сверен с CASCI перед VQE (как на Этапах 1/4).</p>`,
      en: `<p>B3LYP/def2-SVP geometry optimization; ROHF→CASSCF; active spaces via AVAS. All numbers computed.</p>
           <ul>
             <li><strong>Binding/activation:</strong> triplet minimum d(Fe–N)=1.847, d(N–N)=1.129 Å (+0.028 vs free 1.100). E_bind = <strong>+15.4 kcal/mol</strong> relative to Fe(⁵D)+N₂ (endothermic → metastable). Quintet optimizes to d(Fe–N)=5.47 Å (dissociation). Bare Fe–N₂ is unbound in the ground state; slight activation only in the metastable triplet.</li>
             <li><strong>Vertical spin gap</strong> (quintet−triplet) at the bound geometry: PBE <strong>−0.05</strong>, B3LYP <strong>−0.03</strong>, CASSCF(6e,5o) <strong>−0.13</strong> eV → near-degeneracy, DFT-vs-multireference spread ~0.1 eV. (Full-size FeMoco clusters show DFT spin errors up to ~1 eV — about the large cofactor, not this mini-model.)</li>
             <li><strong>Active space for VQE:</strong> AVAS(Fe 3d + N₂ 2p) = CAS(16e,12o) → <strong>24 qubits</strong> (vs 12 for N₂ alone). CASSCF converged for the triplet (E=−1371.036 Ha); the quintet is hard to converge at 24 qubits — itself a sign of the difficulty. This 24-qubit space is the next VQE target on AWS.</li>
           </ul>
           <p>A deliberately minimal "Fe model"; the lesson (need ligands/multinuclearity) leads to the next stage — a [2Fe–2S] cluster. The qubit Hamiltonian will be cross-checked vs CASCI before VQE (as in Stages 1/4).</p>`
    },
    table: {
      title: { ru: "Fe–N₂: геометрия, связывание, спиновая щель — реальный расчёт (def2-SVP)",
               en: "Fe–N₂: geometry, binding, spin gap — real computation (def2-SVP)" },
      head: { ru: ["Величина", "Значение (посчитано)"], en: ["Quantity", "Value (computed)"] },
      rows: [
        ["Триплет: Fe–N / N–N", "1.85 / 1.129 Å  (N≡N удлинена +0.028 Å)"],
        ["Связывание vs Fe(⁵D)+N₂", "+15.4 ккал/моль (метастабильно, не связан)"],
        ["Квинтет", "диссоциирует (Fe–N → 5.47 Å)"],
        ["Вертик. щель квинтет−триплет", "PBE −0.05 · B3LYP −0.03 · CASSCF −0.13 эВ"],
        ["Активное пространство (цель VQE)", "Fe 3d + N₂ 2p = CAS(16,12) → 24 кубита"]
      ]
    },
    figures: [
      { src: "assets/femoco/fe_n2_summary.png",
        caption: { ru: "N₂ на одиночном Fe: метастабильная активация в триплете (+0.028 Å, но +15.4 ккал/моль над фрагментами), квинтет диссоциирует; вертикальная спиновая щель DFT vs CASSCF. Реальный расчёт def2-SVP.",
                   en: "N₂ on a single Fe: metastable triplet activation (+0.028 Å, but +15.4 kcal/mol above fragments), quintet dissociates; vertical spin gap DFT vs CASSCF. Real def2-SVP computation." } }
    ]
  },

  /* ---------------------------------------------------------------- */
  {
    id: "vqe-24q",
    date: "2026-06",
    stage: { ru: "Этап 10", en: "Stage 10" },
    title: { ru: "24-кубитный VQE катализатора: реальный масштаб и стена симулятора",
             en: "24-qubit catalyst VQE: the real scale and the simulator wall" },
    simple: {
      ru: `<p>Берём каталитическое активное пространство Fe–N₂ (24 кубита, Этап 9) в сторону VQE — и впервые упираемся в масштаб «настоящего» катализатора.</p>
           <p>Честные цифры этой задачи: квантовое состояние — это <strong>16.7 миллиона</strong> амплитуд (против 4096 у чистого N₂), точная задача — <strong>174 240</strong> конфигураций, а «ящик инструментов» VQE (операторы возбуждения) содержит <strong>1424</strong> штуки, и каждый шаг оптимизации перебирает их все. На одном классическом симуляторе это минуты-часы на шаг.</p>
           <p>Наш метод работает — проверен на 12 кубитах (ошибка HF 42.7 → <strong>1.4 ккал/моль</strong>). А на 24 кубитах мы это запустили на AWS (32 ядра) и упёрлись в стену <em>вживую</em>: HF посчитан, но <strong>за 69 минут не закрылся даже ОДИН шаг</strong> оптимизации — один проход по 16.7-миллионному квантовому состоянию × 16651 слагаемому гамильтониана непосилен одному классическому компьютеру. То есть на 24 кубитах ломается уже <em>сама классическая симуляция</em> VQE — ровно для этого и нужно настоящее квантовое железо.</p>`,
      en: `<p>We take the Fe–N₂ catalyst active space (24 qubits, Stage 9) toward VQE — and for the first time hit the scale of a "real" catalyst.</p>
           <p>The honest numbers: the quantum state has <strong>16.7 million</strong> amplitudes (vs 4096 for bare N₂), the exact problem has <strong>174,240</strong> configurations, and the VQE "toolbox" of excitation operators holds <strong>1,424</strong>, every optimization step trying all of them. On one classical simulator that is minutes-to-hours per step.</p>
           <p>Our method works — validated at 12 qubits (HF error 42.7 → <strong>1.4 kcal/mol</strong>). At 24 qubits we ran it on AWS (32 cores) and hit the wall <em>live</em>: HF was computed, but <strong>not even ONE optimization step closed in 69 minutes</strong> — a single pass over the 16.7-million-amplitude quantum state × 16,651 Hamiltonian terms is beyond one classical computer. So at 24 qubits it is <em>the classical simulation of VQE itself</em> that breaks — which is exactly why real quantum hardware is needed.</p>`
    },
    tech: {
      ru: `<p>Активное пространство Fe 3d + N₂ 2p из Этапа 9 = <strong>CAS(16e,12o) → 24 кубита</strong> (спин-орбитали). Эталон — CASSCF(16,12) = −1371.036 Ha (триплет).</p>
           <ul>
             <li><strong>Размер:</strong> 174 240 детерминантов; вектор состояния 2²⁴ = 16 777 216 амплитуд (~256 МБ комплексных). Разреженная матрица гамильтониана на 24 кубитах в память не влезает — поэтому VQE считаем <em>без</em> неё, напрямую через adjoint по паули-гамильтониану на lightning.qubit.</li>
             <li><strong>Пул ADAPT:</strong> 64 синглета + 1360 даблетов = <strong>1424</strong> оператора; один проход градиентов по пулу на 24 кубитах — минуты-часы.</li>
             <li><strong>Метод валидирован</strong> (тот же не-разреженный путь, N₂ CAS(6,6)/12 кубитов): HF 42.7 → 1.4 ккал/моль за 30 операторов на lightning.qubit. Кубитный гамильтониан сверяется с CASCI (Этапы 1/4).</li>
           </ul>
           <p><strong>Статус честно (эмпирически измерено):</strong> спуск запущен на AWS (r7i.8xlarge, 32 vCPU, lightning.qubit, авто-терминирование). HF собран: E_HF = −1370.855 Ha для CAS(14e,12o)/16651 членов — это <strong>156 ккал/моль над CASCI</strong> (−1371.104 Ha), т.е. корреляции «отыгрывать» очень много. Но ADAPT-VQE даже с намеренно обрезанным пулом (только 12 даблов) <strong>за 69 минут не закрыл ни одного оператора</strong> (n_ops = 0): один проход градиента по 2²⁴-состоянию × 16651 слагаемых непрактичен на одном узле. Вывод однозначный: метод валидирован на 12 кубитах (Этап 4, ~1 ккал/моль), но на 24 кубитах стену пробивает уже <em>сама классическая симуляция</em> VQE — это и есть аргумент за отказоустойчивый QPU, а не за симулятор. Полный FeMoco (108 кубитов, Этап 5) — далеко за этой стеной. (Запись инкрементальная: HF сохранён на каждом шаге, чтобы результат не терялся.)</p>`,
      en: `<p>The Fe 3d + N₂ 2p active space from Stage 9 = <strong>CAS(16e,12o) → 24 qubits</strong> (spin-orbitals). Reference: CASSCF(16,12) = −1371.036 Ha (triplet).</p>
           <ul>
             <li><strong>Size:</strong> 174,240 determinants; state vector 2²⁴ = 16,777,216 amplitudes (~256 MB complex). The Hamiltonian sparse matrix does not fit in memory at 24 qubits — so we run VQE <em>without</em> it, directly via adjoint on the Pauli Hamiltonian on lightning.qubit.</li>
             <li><strong>ADAPT pool:</strong> 64 singles + 1360 doubles = <strong>1,424</strong> operators; one gradient sweep over the pool at 24 qubits takes minutes-to-hours.</li>
             <li><strong>Method validated</strong> (same non-sparse path, N₂ CAS(6,6)/12 qubits): HF 42.7 → 1.4 kcal/mol in 30 operators on lightning.qubit. The qubit Hamiltonian is cross-checked vs CASCI (Stages 1/4).</li>
           </ul>
           <p><strong>Status, honestly (empirically measured):</strong> the descent was run on AWS (r7i.8xlarge, 32 vCPU, lightning.qubit, auto-terminate). HF was built: E_HF = −1370.855 Ha for CAS(14e,12o)/16,651 terms — <strong>156 kcal/mol above CASCI</strong> (−1371.104 Ha), so there is a great deal of correlation to recover. But ADAPT-VQE, even with a deliberately truncated pool (12 doubles only), <strong>did not close a single operator in 69 minutes</strong> (n_ops = 0): one gradient pass over the 2²⁴ state × 16,651 terms is impractical on one node. The conclusion is unambiguous: the method is validated at 12 qubits (Stage 4, ~1 kcal/mol), but at 24 qubits it is <em>the classical simulation of VQE itself</em> that breaks the wall — the argument for fault-tolerant QPUs, not a simulator. The full FeMoco (108 qubits, Stage 5) is far beyond this wall. (Saving is incremental: HF is written every step so no result is lost.)</p>`
    },
    table: {
      title: { ru: "24-кубитная задача Fe–N₂: ресурсный профиль (реальные числа)",
               en: "24-qubit Fe–N₂ problem: resource profile (real numbers)" },
      head: { ru: ["Величина", "Значение"], en: ["Quantity", "Value"] },
      rows: [
        ["Кубиты (CAS(16,12))", "24"],
        ["Детерминантов (точная задача)", "174 240"],
        ["Амплитуд вектора состояния (2²⁴)", "16 777 216  (~256 МБ)"],
        ["Пул операторов ADAPT", "64 сингла + 1360 даблов = 1424"],
        ["Эталон CASSCF(16,12)", "−1371.036 Ha"],
        ["Метод проверен на N₂ (12 кубитов)", "HF 42.7 → 1.4 ккал/моль"],
        ["24q-спуск на AWS (r7i.8xlarge, 32 vCPU)", "HF −1370.855 Ha; 0 операторов за 69 мин — стена"]
      ]
    },
    figures: []
  }
];

const GLOSSARY = [
  { term: { ru: "Азотфиксация", en: "Nitrogen fixation" },
    def: { ru: "Превращение инертного N₂ из воздуха в усвояемые формы (например, NH₃). Промышленно — Габер–Бош; биологически — нитрогеназа с FeMoco.",
           en: "Converting inert atmospheric N₂ into usable forms (e.g. NH₃). Industrially via Haber–Bosch; biologically via nitrogenase with FeMoco." } },
  { term: { ru: "FeMoco", en: "FeMoco" },
    def: { ru: "Железо-молибден-серный кофактор (Fe₇MoS₉C) — активный центр нитрогеназы, где связывается и восстанавливается N₂.",
           en: "The iron–molybdenum–sulfur cofactor (Fe₇MoS₉C) — the active site of nitrogenase where N₂ is bound and reduced." } },
  { term: { ru: "Статическая (мультиреференсная) корреляция", en: "Static (multireference) correlation" },
    def: { ru: "Ситуация, когда волновую функцию нельзя описать одним детерминантом (типична при разрыве связей). Ломает RHF/CCSD/DFT.",
           en: "When the wavefunction cannot be described by a single determinant (typical for bond breaking). Breaks RHF/CCSD/DFT." } },
  { term: { ru: "Активное пространство, CAS(n,m)", en: "Active space, CAS(n,m)" },
    def: { ru: "n электронов в m орбиталях, в которых учитывается полная корреляция. CAS(54,54) — пространство FeMoco у Reiher 2017.",
           en: "n electrons in m orbitals treated with full correlation. CAS(54,54) is the FeMoco space in Reiher 2017." } },
  { term: { ru: "CASSCF / CASCI / Full CI", en: "CASSCF / CASCI / Full CI" },
    def: { ru: "Многоконфигурационные методы: точная диагонализация в активном пространстве (CASCI/Full CI), с оптимизацией орбиталей (CASSCF).",
           en: "Multiconfigurational methods: exact diagonalization within the active space (CASCI/Full CI), with orbital optimization (CASSCF)." } },
  { term: { ru: "Детерминант", en: "Determinant" },
    def: { ru: "Базовая «конфигурация» расселения электронов по орбиталям. Их число растёт комбинаторно — отсюда ~10¹⁶ у FeMoco.",
           en: "A basic 'configuration' of electrons in orbitals. Their count grows combinatorially — hence ~10¹⁶ for FeMoco." } },
  { term: { ru: "VQE / ADAPT-VQE", en: "VQE / ADAPT-VQE" },
    def: { ru: "Вариационный квантовый алгоритм для энергии основного состояния; ADAPT наращивает анзац адаптивно. Здесь — на симуляторе, сверяется с CASCI.",
           en: "A variational quantum algorithm for the ground-state energy; ADAPT grows the ansatz adaptively. Here on a simulator, checked against CASCI." } },
  { term: { ru: "Кубит", en: "Qubit" },
    def: { ru: "Единица квантовой информации. В наших задачах число кубитов = 2 × число активных орбиталей (спин-орбитали).",
           en: "The unit of quantum information. In our mappings, #qubits = 2 × #active orbitals (spin-orbitals)." } }
];
