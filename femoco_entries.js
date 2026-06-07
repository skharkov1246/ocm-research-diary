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
    id: "scaling",
    date: "2026-06",
    stage: { ru: "Этап 4", en: "Stage 4" },
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
