/* ====================================================================
   OCM research diary — content (RU/EN × simple/technical).
   To add a new stage: append one object to ENTRIES and drop figures in assets/.
   HTML is allowed in the content strings (we author it, so innerHTML is safe).
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
    ru: "Принцип проекта: каждое число — из реального расчёта. Что нельзя проверить — помечаем явно. Ничего не приукрашиваем; неудачи и артефакты показываем.",
    en: "Project rule: every number comes from a real computation. What can't be validated is labelled. Nothing is smoothed; failures and artefacts are shown."
  }
};

const SITE = {
  title:    { ru: "MatterForge diary", en: "MatterForge diary" },
  subtitle: { ru: "Квантовая разработка катализаторов — открытый честный дневник",
              en: "Quantum catalyst design — an open, honest diary" },
  hero: {
    simple: {
      ru: "Мы учим компьютер считать поведение электронов, чтобы понять, можно ли превратить природный газ (метан) сразу в этилен — вещество, из которого делают почти весь пластик. В один шаг это не получается уже 40 лет. Здесь я после каждого этапа честно записываю, что мы посчитали и что узнали.",
      en: "We teach a computer to model how electrons behave, to find out whether natural gas (methane) can be turned in one step into ethylene — the molecule almost all plastics are made from. Nobody has managed this in one step for 40 years. After each stage I honestly log here what we computed and what we learned."
    },
    tech: {
      ru: "Окислительная конденсация метана (OCM) упирается в ~28% потолок выхода C₂ за проход. Мы моделируем лимитирующую стадию (отрыв атома H радикальным кислородом катализатора) методами квантовой химии (CASCI / VQE), строим количественную «карту потолка» селективности и честно проверяем гипотезы о том, как его пробить.",
      en: "Oxidative coupling of methane (OCM) hits a ~28% single-pass C₂-yield ceiling. We model the rate-limiting step (H-atom transfer to the catalyst's radical oxygen) with quantum chemistry (CASCI / VQE), build a quantitative selectivity 'ceiling map', and honestly test hypotheses for breaking it."
    }
  }
};

const ENTRIES = [
  /* ---------------------------------------------------------------- */
  {
    id: "why",
    date: "2026-05",
    stage: { ru: "Зачем", en: "Why" },
    title: { ru: "Зачем всё это: метан, этилен и стена в 28%",
             en: "Why this matters: methane, ethylene, and the 28% wall" },
    simple: {
      ru: `<p>Метан — это природный газ. Этилен — это «кирпичик», из которого делают полиэтилен, то есть почти весь пластик в мире. Сегодня, чтобы получить этилен, газ возят, сжижают, гоняют через много дорогих установок.</p>
           <p><strong>Мечта:</strong> один аппарат, который берёт метан и кислород и сразу выдаёт этилен. Это называется <em>окислительная конденсация метана (OCM)</em>.</p>
           <p><strong>Проблема:</strong> тот же катализатор, который «откусывает» водород у метана, ещё быстрее разрушает уже готовый этилен — сжигает его в углекислый газ. Поэтому за один проход больше ~28% метана в этилен не превратить. В эту стену учёные бьются уже 40 лет.</p>
           <p><strong>Наш подход:</strong> не угадывать катализаторы вслепую, а посчитать на компьютере физику этой стены — и понять, где у неё слабое место.</p>`,
      en: `<p>Methane is natural gas. Ethylene is the building block of polyethylene — almost all the world's plastic. Today, making ethylene means shipping, liquefying and routing gas through many expensive plants.</p>
           <p><strong>The dream:</strong> one reactor that takes methane + oxygen and directly outputs ethylene. This is <em>oxidative coupling of methane (OCM)</em>.</p>
           <p><strong>The problem:</strong> the same catalyst that "bites" a hydrogen off methane destroys the finished ethylene even faster — burning it to CO₂. So a single pass can't convert more than ~28% of methane to ethylene. Scientists have hit this wall for 40 years.</p>
           <p><strong>Our approach:</strong> don't guess catalysts blindly — compute the physics of the wall and find where it is weak.</p>`
    },
    tech: {
      ru: `<p>Реакция мечты: <code>2 CH₄ + ½ O₂ → C₂H₄ + H₂O</code>. Лимитирующая стадия — отрыв атома водорода (HAT) от C–H связи метана радикальным кислородом катализатора (например, центр [Li⁺O⁻] на MgO, классика Лансфорда).</p>
           <p><strong>Почему потолок ~28%:</strong> работает <em>масштабирующее соотношение Бренстеда–Эванса–Полани (BEP)</em>. C–H связи в продуктах C₂ (этан 420, этилен ещё активнее по аллильному/винильному типу) слабее, чем в метане (439 кДж/моль). Любой центр, активный к метану, ещё активнее к продуктам → последовательное переокисление C₂ → COₓ. Это не дефект конкретного катализатора, а свойство механизма.</p>
           <p><strong>Зачем квантовая химия:</strong> переходные состояния с разрывом связей сильно мультиреференсны — однодетерминантные методы (обычный DFT/HF) систематически врут. Нужны многоконфигурационные методы (CASCI/CASSCF), а в перспективе — VQE на квантовом компьютере для активных пространств, недоступных классике (CAS ≳ 20 орбиталей).</p>`,
      en: `<p>Dream reaction: <code>2 CH₄ + ½ O₂ → C₂H₄ + H₂O</code>. The rate-limiting step is hydrogen-atom transfer (HAT) from a methane C–H bond to the catalyst's radical oxygen (e.g. the [Li⁺O⁻] centre on MgO, the classic Lunsford site).</p>
           <p><strong>Why the ~28% ceiling:</strong> the <em>Brønsted–Evans–Polanyi (BEP) scaling relation</em>. The C–H bonds of the C₂ products (ethane 420; ethylene even more reactive via allylic/vinylic paths) are weaker than methane's (439 kJ/mol). Any site active toward methane is more active toward the products → sequential over-oxidation C₂ → COₓ. This is a property of the mechanism, not a flaw of one catalyst.</p>
           <p><strong>Why quantum chemistry:</strong> bond-breaking transition states are strongly multireference — single-determinant methods (ordinary DFT/HF) are systematically wrong. We need multiconfigurational methods (CASCI/CASSCF), and eventually VQE on a quantum computer for active spaces beyond classical reach (CAS ≳ 20 orbitals).</p>`
    },
    figures: []
  },

  /* ---------------------------------------------------------------- */
  {
    id: "phase1",
    date: "2026-05",
    stage: { ru: "Этап 1", en: "Stage 1" },
    title: { ru: "Честный калькулятор: собрали и проверили движок",
             en: "An honest calculator: building and validating the engine" },
    simple: {
      ru: `<p>Сначала мы построили «честный калькулятор» — программу, которая считает энергию молекул на каждом шаге реакции и рисует «горку», которую нужно преодолеть (это и есть энергетический барьер).</p>
           <p>Главное правило: <strong>никакого вранья</strong>. Если квантовый метод не сошёлся или число выглядит нефизично — программа громко об этом говорит, а не подгоняет результат.</p>
           <p>Мы проверили калькулятор на простых молекулах (H₂, LiH), где ответ известен точно: совпадение — до 13-го знака. Значит, фундамент надёжный.</p>
           <p>На графике слева — как «квантовый» алгоритм (VQE) шаг за шагом приближается к точному ответу. Справа — энергетическая «горка» реакции.</p>`,
      en: `<p>First we built an "honest calculator" — a program that computes the energy of the molecules at each step of the reaction and draws the "hill" that must be climbed (the energy barrier).</p>
           <p>The core rule: <strong>no fibbing</strong>. If a quantum method fails to converge, or a number looks unphysical, the program says so loudly instead of fudging the result.</p>
           <p>We checked the calculator on simple molecules (H₂, LiH) with exactly known answers: agreement to the 13th digit. So the foundation is solid.</p>
           <p>The figure on the left shows how the "quantum" algorithm (VQE) steps closer to the exact answer. The reaction's energy "hill" is on the right.</p>`
    },
    tech: {
      ru: `<p>Собрали воспроизводимый пайплайн на <strong>PennyLane + PySCF</strong>: кластер Li/MgO, линейная синхронная транзитная (LST) траектория HAT-стадии, активное пространство <strong>CAS(7e,6o) → 12 кубитов</strong> (Jordan–Wigner), и сравнение VQE с точным CASCI.</p>
           <p><strong>Валидация рецепта</strong> «интегралы → кубитный гамильтониан»: на H₂ и LiH основное состояние совпадает с CASCI до <strong>~1e-13 ккал/моль</strong> (правильная конвенция ERI, транспонирование (0,2,3,1)). CI-гейт стоит перед smoke-тестом.</p>
           <p><strong>Честный результат VQE</strong> на сильно-коррелированном CAS(7,6): UCCSD сходится к своему вариационному минимуму ~1.3–1.5 ккал/моль <em>выше</em> CASCI (singles+doubles не достают до FCI здесь); ADAPT-VQE за 5-мин бюджет дошёл до 7.8 ккал/моль в реактанте. Критерий «<1 ккал/моль» выполнен чисто на контрольном CAS(4,4) (0.21 ккал/моль). Это сообщается честно, без сглаживания.</p>
           <p><strong>Производный барьер C–H активации</strong>: E(TS)−E(реактант) = <strong>74.35 ккал/моль</strong> (CASCI; приближённый TS при λ=0.62).</p>`,
      en: `<p>We built a reproducible pipeline on <strong>PennyLane + PySCF</strong>: a Li/MgO cluster, a linear-synchronous-transit (LST) path for the HAT step, an active space of <strong>CAS(7e,6o) → 12 qubits</strong> (Jordan–Wigner), and VQE benchmarked against exact CASCI.</p>
           <p><strong>Recipe validation</strong> "integrals → qubit Hamiltonian": on H₂ and LiH the ground state matches CASCI to <strong>~1e-13 kcal/mol</strong> (correct ERI convention, transpose (0,2,3,1)). This CI gate runs before the smoke test.</p>
           <p><strong>Honest VQE result</strong> on the strongly-correlated CAS(7,6): UCCSD converges to its variational minimum ~1.3–1.5 kcal/mol <em>above</em> CASCI (singles+doubles can't reach FCI here); ADAPT-VQE reached 7.8 kcal/mol at the reactant within a 5-min budget. The "&lt;1 kcal/mol" criterion is met cleanly on the CAS(4,4) control (0.21 kcal/mol). Reported honestly, not smoothed.</p>
           <p><strong>Derived C–H activation barrier</strong>: E(TS)−E(reactant) = <strong>74.35 kcal/mol</strong> (CASCI; approximate TS at λ=0.62).</p>`
    },
    figures: [
      { src: "assets/vqe_convergence.png",
        caption: { ru: "Сходимость VQE к точному (CASCI) ответу по ходу оптимизации.",
                   en: "VQE convergence toward the exact (CASCI) answer during optimization." } },
      { src: "assets/reaction_profile.png",
        caption: { ru: "Энергетический профиль HAT-стадии: барьер ≈ 74 ккал/моль (CASCI).",
                   en: "Energy profile of the HAT step: barrier ≈ 74 kcal/mol (CASCI)." } },
      { src: "assets/scaling.png",
        caption: { ru: "Рост классической стоимости с размером активного пространства (НЕ заявление о квантовом превосходстве).",
                   en: "Classical cost vs. active-space size (NOT a quantum-advantage claim)." } }
    ]
  },

  /* ---------------------------------------------------------------- */
  {
    id: "engine",
    date: "2026-06-04",
    stage: { ru: "Этап 2", en: "Stage 2" },
    title: { ru: "Движок открытий и карта потолка: насколько высоко 70%?",
             en: "The discovery engine & the ceiling map: how far is 70%?" },
    simple: {
      ru: `<p>Чтобы не перебирать катализаторы наугад, мы придумали <strong>один числовой «индикатор селективности»</strong>. Он сравнивает: насколько легко центр атакует метан против того, насколько легко он атакует продукт.</p>
           <ul><li><strong>Минус</strong> = центр быстрее портит продукт → плохо (как раз наша стена).</li>
               <li><strong>Плюс</strong> = центр бережёт продукт → вот это бы пробило стену.</li></ul>
           <p>Проверка честности: на знаменитом катализаторе Li/MgO индикатор вышел <strong>отрицательным</strong> — ровно как в реальности (он действительно переокисляет). Значит, движок ловит правду. Тот же расчёт мы повторили в облаке (AWS) — получили то же число.</p>
           <p>Дальше мы построили <strong>карту потолка</strong>: какой индикатор нужен для какого выхода. Вердикт: для мечты в <strong>70%</strong> нужен индикатор около <strong>+4</strong>. А все «электронные» подкрутки катализатора дают максимум −0.6 — то есть остаются у стены. Электроникой стену не пробить.</p>`,
      en: `<p>Instead of trying catalysts at random, we invented <strong>one numeric "selectivity indicator"</strong>. It compares how easily a site attacks methane versus how easily it attacks the product.</p>
           <ul><li><strong>Negative</strong> = the site ruins the product faster → bad (exactly our wall).</li>
               <li><strong>Positive</strong> = the site spares the product → this would break the wall.</li></ul>
           <p>Honesty check: on the famous Li/MgO catalyst the indicator came out <strong>negative</strong> — exactly matching reality (it really does over-oxidize). So the engine catches the truth. We repeated the same computation in the cloud (AWS) and got the same number.</p>
           <p>Then we built a <strong>ceiling map</strong>: which indicator value is needed for which yield. Verdict: the <strong>70%</strong> dream needs an indicator around <strong>+4</strong>. But all "electronic" tweaks of the catalyst give at most −0.6 — still stuck at the wall. Electronics alone can't break it.</p>`
    },
    tech: {
      ru: `<p><strong>Дескриптор:</strong> <code>ΔΔG‡ = барьер(C₂H₆ C–H) − барьер(CH₄ C–H)</code> на одном центре, одинаковым жёстко-каркасным HAT-протоколом для обоих субстратов → <em>разность</em> устойчива к систематической ошибке. Движок: ROHF + frontier-CAS CASCI / STO-3G (UKS-DFT на открытой оболочке оксида не сходился; CASCI — сходится).</p>
           <p><strong>Гейт честности:</strong> Li/MgO O⁻ — известный переокислитель → ждём ΔΔG‡ &lt; 0. Получили <strong>−1.09 ккал/моль</strong>, воспроизведено на AWS (узел вернул −1.088). Пайплайн корректен и портируем.</p>
           <p><strong>Скрин по электронике/геометрии</strong> (допанты, напряжение Li–O, сближение C···O): лучшее ΔΔG‡ ≈ <strong>−0.60</strong> — остаётся BEP-связанным. Два экстремальных геометрии дали SCF-артефакты — помечены и исключены.</p>
           <p><strong>Карта потолка:</strong> ΔΔG‡ → выход через модель последовательной реакции CH₄→C₂→COₓ, Y_max = κ^(κ/(1−κ)), κ = C·exp(−ΔΔG‡/RT), калибровка на Li/MgO ≈ 25%. Устойчиво по T и калибровке: для 70% нужно <strong>ΔΔG‡ ≈ +4 (диапазон +3.1…+5.0)</strong>.</p>`,
      en: `<p><strong>Descriptor:</strong> <code>ΔΔG‡ = barrier(C₂H₆ C–H) − barrier(CH₄ C–H)</code> on one site, with the same rigid-frame HAT protocol for both substrates → the <em>difference</em> is robust to systematic error. Engine: ROHF + frontier-CAS CASCI / STO-3G (UKS-DFT stalled on the open-shell oxide; CASCI converges).</p>
           <p><strong>Honesty gate:</strong> Li/MgO O⁻ is a known over-oxidizer → expect ΔΔG‡ &lt; 0. We get <strong>−1.09 kcal/mol</strong>, reproduced on AWS (node returned −1.088). The pipeline is correct and portable.</p>
           <p><strong>Electronic/geometric screen</strong> (dopants, Li–O strain, C···O approach): best ΔΔG‡ ≈ <strong>−0.60</strong> — stays BEP-bound. Two extreme geometries gave SCF artefacts — flagged and excluded.</p>
           <p><strong>Ceiling map:</strong> ΔΔG‡ → yield via a consecutive-reaction model CH₄→C₂→COₓ, Y_max = κ^(κ/(1−κ)), κ = C·exp(−ΔΔG‡/RT), calibrated to Li/MgO ≈ 25%. Robust across T and calibration: 70% needs <strong>ΔΔG‡ ≈ +4 (range +3.1…+5.0)</strong>.</p>`
    },
    figures: [
      { src: "assets/ceiling_map.png",
        caption: { ru: "Карта потолка: какой ΔΔG‡ нужен для какого выхода C₂. Зелёная зона (ΔΔG‡>0) пуста — ни один электронный центр туда не дотягивает.",
                   en: "Ceiling map: which ΔΔG‡ each C₂ yield requires. The green zone (ΔΔG‡>0) is empty — no electronic site reaches it." } }
    ],
    table: {
      title: { ru: "Сколько нужно для каждого выхода", en: "What each yield demands" },
      head: { ru: ["Цель (выход C₂)", "Нужный ΔΔG‡"], en: ["Target (C₂ yield)", "Required ΔΔG‡"] },
      rows: [
        ["28% (текущий потолок / current ceiling)", "−0.71"],
        ["40%", "+0.67"],
        ["50%", "+1.73"],
        ["70%", "+3.98 ккал/моль"]
      ]
    }
  },

  /* ---------------------------------------------------------------- */
  {
    id: "stresstest",
    date: "2026-06-05",
    stage: { ru: "Этап 3", en: "Stage 3" },
    title: { ru: "Честный стресс-тест: как движок поймал собственный оптимизм",
             en: "The honest stress-test: how the engine caught its own optimism" },
    simple: {
      ru: `<p>У нас была красивая идея «двери из стены»: посадить активный кислород на дно крошечной пещерки (поры). Метан маленький — заходит и реагирует; продукт побольше — не может вернуться, чтобы его сожгли. Это «молекулярное сито».</p>
           <p>Первый прикидочный расчёт (на этане) дал огромный запас: +13.6 при нужных +4. «Достаточно!»</p>
           <p><strong>Но мы проверили честно</strong> — и идея почти рассыпалась:</p>
           <ul>
             <li>защищать надо не этан, а <strong>этилен</strong> — а он по размеру <em>почти как метан</em>;</li>
             <li>молекула может повернуться поудобнее и проскользнуть в пору.</li>
           </ul>
           <p>Итог: пора, в которую входит метан, <strong>впускает и этилен тоже</strong>. Разделить их по размеру не выходит. Лучшее, что даёт чистое сито — около <strong>42%</strong>, а не 70%.</p>
           <p>Это и есть смысл проекта: <strong>движок поймал и исправил собственную ошибку</strong>. А ещё мы перепроверили главный индикатор на всё более точной математике — вывод не перевернулся. Фундамент крепкий, а красивую идею пришлось честно урезать.</p>`,
      en: `<p>We had a pretty "door in the wall" idea: put the reactive oxygen at the bottom of a tiny cave (a pore). Methane is small — it enters and reacts; the bigger product can't get back in to be burned. This is "molecular sieving".</p>
           <p>A first rough estimate (on ethane) gave a huge margin: +13.6 against the +4 needed. "Sufficient!"</p>
           <p><strong>But we checked it honestly</strong> — and the idea nearly fell apart:</p>
           <ul>
             <li>the molecule to protect is not ethane but <strong>ethylene</strong> — which is <em>almost the same size as methane</em>;</li>
             <li>a molecule can rotate into a comfier pose and slip into the pore.</li>
           </ul>
           <p>Result: a pore that admits methane <strong>admits ethylene too</strong>. You can't separate them by size. The best a pure sieve gives is about <strong>42%</strong>, not 70%.</p>
           <p>This is the whole point of the project: <strong>the engine caught and corrected its own error</strong>. We also re-checked the key indicator with progressively more accurate math — the conclusion didn't flip. The foundation is solid; the pretty idea had to be honestly trimmed.</p>`
    },
    tech: {
      ru: `<p><strong>Гипотеза:</strong> конфайнмент-сайт ломает BEP через форм-селективность, а не через силу связи. Стерический рычаг считаем закрытооболочечной RHF/STO-3G в He-клетке (чистый паули-потенциал, без артефактов открытой оболочки).</p>
           <p><strong>Первая оценка (sieving.py):</strong> CH₄ vs <em>этан</em>, одна ориентация → ΔΔE_sieve = +13.6 ккал/моль при поре 3.2 Å. Выглядело «достаточно» (≫ +4).</p>
           <p><strong>Стресс-тест (confine_products.py)</strong> исправил две ошибки: (1) защищать надо <em>этилен</em> (кин. диаметр ~3.9 Å ≈ метан ~3.8), а не этан (~4.4); (2) энергию конфайнмента надо <strong>минимизировать по ориентации</strong> (молекула входит в самой удобной позе). Результат: при 3.2 Å метан проходит (1.9), этилен <em>тоже</em> (3.9) → дифференциал всего <strong>+1.9</strong>. Этилен ≈3× легче впустить, чем этан (отношение ≈0.32). Чистого окна «метан-внутрь / этилен-наружу» нет: при 2.9 Å этилен наконец режется (+6.2), но метан уже сам душится (+7.3).</p>
           <p><strong>На карте потолка (confined_map.py):</strong> ΔΔG‡_eff = −1.09 + η·ΔΔE_access. С честным +1.9 даже при η=1 выход <strong>~42%, не 70%</strong>.</p>
           <p><strong>Устойчивость фундамента (descriptor_basis.py):</strong> знак ΔΔG‡ по базисам sto-3g/3-21g/6-31g/6-31g* = −1.09 / −3.20 / −3.67 / −2.00 — <strong>везде &lt; 0</strong>, все SCF сошлись. Неселективность — не артефакт минимального базиса.</p>
           <p><strong>Вторая дверь — неокислительный путь</strong> (2CH₄→C₂H₄+2H₂): без O₂ нет COₓ-потолка, но равновесно-ограничен; ~22%/41% при 1090 °C (1 бар / отвод H₂) — режим Bao. Для 70% нужно &gt;1200 °C.</p>
           <p><strong>Честный итог:</strong> 70% в один шаг не даёт <em>ни один</em> просчитанный механизм. Реальные двери требуют большего, чем «подбор сайта»: (A) различение по <em>форме/ориентации</em> (не по размеру) — проверяется только на реальной релаксированной поре (AWS-уровень); (B) неокислительный + отвод H₂.</p>`,
      en: `<p><strong>Hypothesis:</strong> a confined site breaks BEP via shape-selectivity, not bond strength. We quantify the steric lever with closed-shell RHF/STO-3G in a He cage (pure Pauli potential, no open-shell artefacts).</p>
           <p><strong>First estimate (sieving.py):</strong> CH₄ vs <em>ethane</em>, one fixed orientation → ΔΔE_sieve = +13.6 kcal/mol at a 3.2 Å pore. Looked "sufficient" (≫ +4).</p>
           <p><strong>Stress-test (confine_products.py)</strong> fixed two errors: (1) the species to protect is <em>ethylene</em> (kinetic diameter ~3.9 Å ≈ methane ~3.8), not ethane (~4.4); (2) the confinement energy must be <strong>minimized over orientations</strong> (a molecule enters in its best-fit pose). Result: at 3.2 Å methane is admitted (1.9) and ethylene is admitted <em>too</em> (3.9) → differential only <strong>+1.9</strong>. Ethylene is ~3× easier to admit than ethane (ratio ≈0.32). No clean "methane-in / ethylene-out" window: at 2.9 Å ethylene is finally penalized (+6.2) but methane is itself throttled (+7.3).</p>
           <p><strong>On the ceiling map (confined_map.py):</strong> ΔΔG‡_eff = −1.09 + η·ΔΔE_access. With the honest +1.9, even at η=1 the yield is <strong>~42%, not 70%</strong>.</p>
           <p><strong>Foundation robustness (descriptor_basis.py):</strong> the ΔΔG‡ sign across bases sto-3g/3-21g/6-31g/6-31g* = −1.09 / −3.20 / −3.67 / −2.00 — <strong>negative everywhere</strong>, all SCF converged. Non-selectivity is not a minimal-basis artefact.</p>
           <p><strong>Second door — the non-oxidative lane</strong> (2CH₄→C₂H₄+2H₂): no O₂ means no COₓ ceiling, but it's equilibrium-limited; ~22%/41% at 1090 °C (1 bar / H₂ removal) — the Bao regime. 70% needs &gt;1200 °C.</p>
           <p><strong>Honest bottom line:</strong> <em>no</em> screened mechanism reaches 70% in one step. The real doors need more than site-tuning: (A) discrimination by <em>shape/orientation</em> (not size) — only testable on a real relaxed pore (AWS-scale); (B) non-oxidative + H₂ removal.</p>`
    },
    figures: [
      { src: "assets/confine_products.png",
        caption: { ru: "Стресс-тест: кривые CH₄ и этилена почти совпадают — по размеру их не разделить.",
                   en: "Stress-test: the CH₄ and ethylene curves nearly overlap — they can't be separated by size." } },
      { src: "assets/confined_map.png",
        caption: { ru: "Где оказывается конфайнмент-сайт на карте: при поре, пропускающей метан, потолок ~42%.",
                   en: "Where the confined site lands on the map: at the methane-admitting pore, the cap is ~42%." } },
      { src: "assets/phase3_summary.png",
        caption: { ru: "Итоговая двухпанельная картина: ни электроника, ни чистое сито не дают 70%.",
                   en: "The two-panel summary: neither electronics nor pure sieving reach 70%." } },
      { src: "assets/nonox_thermo.png",
        caption: { ru: "Вторая дверь: неокислительный путь, равновесно-ограниченный, ~41% при 1090 °C с отводом H₂.",
                   en: "The second door: the non-oxidative lane, equilibrium-limited, ~41% at 1090 °C with H₂ removal." } }
    ],
    table: {
      title: { ru: "Стресс-тест сита: энергия входа в пору (ккал/моль, минимум по ориентации)",
               en: "Sieve stress-test: pore-entry energy (kcal/mol, orientation-minimized)" },
      head: { ru: ["Пора R (Å)", "CH₄ (впустить)", "C₂H₄ этилен (отсечь)", "разница"],
              en: ["Pore R (Å)", "CH₄ (admit)", "C₂H₄ ethylene (exclude)", "differential"] },
      rows: [
        ["2.9", "7.3", "13.6", "+6.2 (но метан тоже душится / methane throttled too)"],
        ["3.2", "1.9", "3.9", "+1.9 (этилен протекает / ethylene leaks in)"],
        ["3.6", "0.3", "0.6", "+0.3"]
      ]
    }
  },

  /* ---------------------------------------------------------------- */
  {
    id: "thissite",
    date: "2026-06-05",
    stage: { ru: "Этап 4", en: "Stage 4" },
    title: { ru: "Этот дневник: чтобы учиться вместе с экспериментами",
             en: "This diary: learning alongside the experiments" },
    simple: {
      ru: `<p>Сегодня мы сделали этот сайт. Идея простая: после каждого этапа я записываю сюда, что мы посчитали, на двух языках и двух уровнях сложности — чтобы ты учился параллельно тому, что мы делаем.</p>
           <p>Переключатели вверху меняют <strong>язык</strong> (RU/EN) и <strong>уровень</strong> (простой/технический). Внизу — словарик терминов.</p>
           <p>Каждый новый этап исследования будет появляться здесь новой записью.</p>`,
      en: `<p>Today we built this site. The idea is simple: after each stage I log here what we computed, in two languages and two difficulty levels — so you learn in parallel with the work.</p>
           <p>The toggles at the top switch <strong>language</strong> (RU/EN) and <strong>level</strong> (simple/technical). A glossary is at the bottom.</p>
           <p>Every new research stage will appear here as a new entry.</p>`
    },
    tech: {
      ru: `<p>Статический сайт (HTML/CSS/JS, без сборки), задеплоен на <strong>GitHub Pages</strong>. Контент — в <code>entries.js</code> как массив записей с полями simple/tech × ru/en; рендер пересобирает страницу при переключении состояния (язык/уровень сохраняются в localStorage).</p>
           <p>Добавить этап = дописать один объект в <code>ENTRIES</code> и положить фигуры в <code>assets/</code>, затем <code>git push</code> — Pages пересоберётся сам. Так дневник растёт вместе с исследованием.</p>`,
      en: `<p>A static site (HTML/CSS/JS, no build step) deployed on <strong>GitHub Pages</strong>. Content lives in <code>entries.js</code> as an array of records with simple/tech × ru/en fields; the renderer rebuilds the page on state change (language/level persist in localStorage).</p>
           <p>Adding a stage = append one object to <code>ENTRIES</code>, drop figures in <code>assets/</code>, then <code>git push</code> — Pages rebuilds itself. The diary grows with the research.</p>`
    },
    figures: []
  },

  /* ---------------------------------------------------------------- */
  {
    id: "confined",
    date: "2026-06-06",
    stage: { ru: "Этап 5", en: "Stage 5" },
    title: { ru: "Конфайнмент на AWS: красивый сигнал, который не выжил",
             en: "Confinement on AWS: a beautiful signal that didn't survive" },
    simple: {
      ru: `<p>Была красивая гипотеза: посадить активный кислород на дно крошечной «пещерки» из
           атомов (поры) и так зажать молекулы, чтобы селективность перевернулась — то, чего не
           смогли ни электроника, ни сито по размеру. Дешёвый расчёт радостно сказал: <strong>ДА,
           +6.3</strong> — выше порога 70%! 🎉</p>
           <p>Но прежде чем поверить, мы сделали <strong>контроль</strong>: отодвинули стенки
           «пещерки» далеко. Логика простая — если эффект от <em>зажатия</em>, то с далёкими
           стенками он должен исчезнуть. И на нормальной математике (на облачном компьютере AWS)
           контроль <strong>провалился</strong>: с далёкими стенками эффект стал ещё
           <strong>больше</strong>. Значит, дело не в зажатии — просто «добавили атомы → число
           поменялось». Ложная тревога.</p>
           <p><strong>Прорыва нет.</strong> Но главное — мы <strong>поймали обманку до того, как в
           неё поверили</strong>. Вот так и работает честная наука: красивый результат сначала
           пытаются убить контролем, и настоящим считается только выживший. Карта потолка
           (Этапы 1–3) остаётся в силе: 70% за один шаг пока не берёт ни один проверенный
           механизм. (Весь расчёт — на AWS, ~$2, инстанс выключен.)</p>`,
      en: `<p>A pretty hypothesis: put the reactive oxygen at the bottom of a tiny atomic "cave"
           (a pore) and squeeze the molecules so selectivity flips — something neither electronics
           nor size-sieving could do. The cheap calculation happily said: <strong>YES, +6.3</strong>
           — above the 70% threshold! 🎉</p>
           <p>But before believing it, we ran a <strong>control</strong>: move the cave walls far
           away. Simple logic — if the effect is from <em>squeezing</em>, far walls should kill it.
           And with proper math (on an AWS cloud computer) the control <strong>failed</strong>: with
           the walls far away the effect got even <strong>bigger</strong>. So it wasn't squeezing —
           just "add atoms → the number changes." A false alarm.</p>
           <p><strong>No breakthrough.</strong> But the point is we <strong>caught the fake-out
           before believing it</strong>. That's how honest science works: you first try to kill a
           pretty result with a control, and only the survivor counts as real. The ceiling map
           (Stages 1–3) stands: 70% in one step is still beaten by no tested mechanism. (All compute
           on AWS, ~$2, instance shut down.)</p>`
    },
    tech: {
      ru: `<p><strong>Канал и движок:</strong> считали на ВАЛИДИРОВАННОМ канале (этан-HAT,
           воспроизводит −1.09) валидированным движком ROHF + frontier-CAS CASCI (жёсткий скан).
           Конфайнмент — <em>связанная</em> MgO-полость (Mg²⁺/O²⁻ на позициях решётки, целыми
           нейтральными блоками → закрытая оболочка), без плавающего He, который ломал прошлый заход.</p>
           <p><strong>Контроль:</strong> те же атомы воротника, отодвинутые в 2× (присутствуют, но не
           зажимают). Если эффект — конфайнмент, дальний воротник НЕ должен воспроизводить сдвиг.</p>
           <p><strong>STO-3G:</strong> bare −1.09 → ближний +6.34, дальний −2.08 (контроль ✅ —
           выглядело реально). <strong>def2-SVP (AWS):</strong> bare −0.01, ближний +8.70 (у CH₄
           барьер схлопнулся в 0.00 — патология), дальний <strong>+12.91 &gt; ближнего</strong> →
           <strong>контроль ❌</strong>: сдвиг — аддитивный артефакт добавленных атомов, не зажатие.</p>
           <p><strong>Доп. флаг:</strong> голый дескриптор просел −1.09 (STO-3G) → −0.01 (def2-SVP) —
           авто-выбор активного пространства не устойчив по базису. Скрипт сам пометил «NOT TRUSTED».</p>
           <p><strong>Что нужно для настоящего теста:</strong> фиксированное/валидированное активное
           пространство (AVAS), периодическая/embedding-модель реальной поры, и правильный канал для
           этилена (π-окисление, не HAT). Это отдельный большой расчёт.</p>`,
      en: `<p><strong>Channel & engine:</strong> computed on the VALIDATED channel (ethane HAT,
           reproduces −1.09) with the validated ROHF + frontier-CAS CASCI rigid-scan engine.
           Confinement = a <em>bonded</em> MgO cavity (Mg²⁺/O²⁻ at lattice positions, whole neutral
           units → closed shell), no floating He (which broke the previous attempt).</p>
           <p><strong>Control:</strong> the same collar atoms moved 2× out (present but not confining).
           If the effect is confinement, the far collar must NOT reproduce the shift.</p>
           <p><strong>STO-3G:</strong> bare −1.09 → near +6.34, far −2.08 (control ✅ — looked real).
           <strong>def2-SVP (AWS):</strong> bare −0.01, near +8.70 (CH₄ barrier collapsed to 0.00 —
           pathological), far <strong>+12.91 &gt; near</strong> → <strong>control ❌</strong>: the
           shift is an additive artefact of adding atoms, not squeezing.</p>
           <p><strong>Extra flag:</strong> the bare descriptor drifted −1.09 (STO-3G) → −0.01
           (def2-SVP) — the auto active-space selection is not basis-robust. The script itself
           flagged "NOT TRUSTED".</p>
           <p><strong>What a real test needs:</strong> a fixed/validated active space (AVAS), a
           periodic/embedded model of an actual pore, and the right channel for ethylene
           (π-oxidation, not HAT). That is a separate, larger calculation.</p>`
    },
    table: {
      title: { ru: "Дескриптор ΔΔG‡ по базису и контролю (ккал/моль)",
               en: "Descriptor ΔΔG‡ across basis & control (kcal/mol)" },
      head: { ru: ["", "голый сайт", "ближний воротник", "дальний контроль"],
              en: ["", "bare site", "near collar", "far control"] },
      rows: [
        ["STO-3G", "−1.09", "+6.34", "−2.08  ✅ control ok"],
        ["def2-SVP", "−0.01", "+8.70", "+12.91  ❌ control fails"]
      ]
    },
    figures: [
      { src: "assets/cavity_selfcorrection.png",
        caption: { ru: "Самокоррекция: STO-3G дал заманчивый +6.3, но дальний контроль на def2-SVP даёт ещё больше (+12.9) — значит это артефакт добавленных атомов, а не конфайнмент.",
                   en: "Self-correction: STO-3G gave a tempting +6.3, but at def2-SVP the far control gives even more (+12.9) — so it is an additive artefact, not confinement." } }
    ]
  },

  /* ---------------------------------------------------------------- */
  {
    id: "routes",
    date: "2026-06-06",
    stage: { ru: "Этап 6", en: "Stage 6" },
    title: { ru: "Третий путь: перебор маршрутов — где революционная магистраль?",
             en: "The third path: screening routes — where is the revolutionary trunk?" },
    simple: {
      ru: `<p>Вместо того чтобы чинить кислородный (O₂) катализатор, мы зашли с совсем другой
           стороны: <strong>перебрали много маршрутов</strong> «метан → этилен» с разными
           «заменителями кислорода» и при разных температурах — хлор, бром, сера, CO₂, вообще без
           окислителя. Зачем? Стена O₂ в том, что кислород <em>сжигает</em> готовый продукт в CO₂.
           А если взять «мягкого» партнёра, который не сжигает?</p>
           <p><strong>Результат:</strong> с <strong>хлором/бромом</strong> метан соединяется почти
           полностью <em>уже при 400 °C</em>, и продукт НЕ горит (нет CO₂) — потому что реакция идёт
           через устойчивый промежуток (CH₃Cl), который сшивают отдельно. <strong>Сера</strong>
           (мягкий окислитель) тоже обходит «горение». А обычный O₂ даёт «100% конверсии» только на
           бумаге — на деле горит в CO₂ и упирается в свои ~28%.</p>
           <p><strong>Магистраль:</strong> революция не в «лучшем O₂-катализаторе», а в <strong>смене
           химии на медиатор, который не жжёт продукт</strong> (галогенный цикл / серный мягкий путь).
           Подвох честно: хлор/серу надо <em>возвращать в цикл</em>, и хлор склонен «пережёвывать»
           свой продукт (бром чище). Но стены COx здесь <strong>механически нет</strong> — это и есть
           неординарный путь.</p>`,
      en: `<p>Instead of fixing the oxygen (O₂) catalyst, we came at it from a completely different
           angle: <strong>screened many routes</strong> for methane → ethylene with different "oxygen
           substitutes" and temperatures — chlorine, bromine, sulfur, CO₂, no oxidant at all. Why? The
           O₂ wall is that oxygen <em>burns</em> the finished product to CO₂. What if we use a "soft"
           partner that doesn't burn it?</p>
           <p><strong>Result:</strong> with <strong>chlorine/bromine</strong> methane couples almost
           completely <em>already at 400 °C</em>, and the product is NOT burned (no CO₂) — because it
           goes through a stable intermediate (CH₃Cl) that is coupled separately. <strong>Sulfur</strong>
           (a soft oxidant) also escapes the "burning". Plain O₂ shows "100% conversion" only on paper —
           in reality it burns to CO₂ and hits its ~28%.</p>
           <p><strong>The trunk:</strong> the revolution is not a "better O₂ catalyst" but
           <strong>switching the chemistry to a mediator that doesn't burn the product</strong>
           (a halogen loop / soft sulfur route). The honest catch: chlorine/sulfur must be
           <em>recycled</em>, and chlorine tends to over-react with its own product (bromine is
           cleaner). But the COx wall is <strong>mechanically absent</strong> here — that's the
           non-ordinary path.</p>`
    },
    tech: {
      ru: `<p><strong>Термодинамическая карта</strong> (NIST ΔHf/S°, равновесный экстент) по 12
           маршрутам × T (400–1400 °C, 1 бар): трунк = <em>мягкие/медиаторные</em> окислители.
           Cl₂/Br₂ → C₂H₆: <strong>~100% / 93% уже при 400 °C, без COx</strong>; S₂ → C₂H₄ ~61% при
           высокой T, без COx; неокислительный/MDA — равновесно-ограничены (~51–61% при 1400 °C); O₂
           термодинамически ~100%, но <strong>COx-стена</strong> режет выход до ~28% (конверсия не
           лимит — лимит селективность).</p>
           <p><strong>Квантовая проверка (газофазный HAT, UKS/PBE0):</strong> дескриптор
           ΔΔG‡ = барьер(X•+C₂H₆) − барьер(X•+CH₄) для X = OH/Cl/Br/SH — <strong>у ВСЕХ &lt; 0</strong>
           (ни один радикал-абстрактор не инвертирует селективность). Это Эванс–Поляни: разница BDE
           C–H субстратов не зависит от абстрактора. <strong>Вывод: преимущество трунка —
           механистическое</strong> (устойчивый CH₃X развязывает активацию и переокисление; мягкая
           термодинамика серы), а не «магический катализатор». Переокисление продукта: Cl хуже Br
           (бромный путь чище). <em>(def2-SVP подтверждение знаков считается на AWS.)</em></p>
           <p><strong>Честные оговорки:</strong> термодинамика = <em>способность</em>, не выход
           (кинетика/коксование/разделение решают остальное); абсолюты жёсткого скана грубые (знаки
           надёжны). Узкое место смещается к <em>петле медиатора</em> (рецикл Cl/HX, управление H₂S)
           и селективности сшивки CH₃X → C₂ — это следующий настоящий AWS-расчёт.</p>`,
      en: `<p><strong>Thermodynamic map</strong> (NIST ΔHf/S°, equilibrium extent) over 12 routes × T
           (400–1400 °C, 1 bar): the trunk is <em>soft/mediated</em> oxidants. Cl₂/Br₂ → C₂H₆:
           <strong>~100% / 93% already at 400 °C, no COx</strong>; S₂ → C₂H₄ ~61% at high T, no COx;
           non-oxidative/MDA are equilibrium-limited (~51–61% at 1400 °C); O₂ is thermodynamically
           ~100% but the <strong>COx wall</strong> caps yield at ~28% (conversion isn't the limit —
           selectivity is).</p>
           <p><strong>Quantum check (gas-phase HAT, UKS/PBE0):</strong> the descriptor
           ΔΔG‡ = barrier(X•+C₂H₆) − barrier(X•+CH₄) for X = OH/Cl/Br/SH is <strong>&lt; 0 for ALL</strong>
           (no radical abstractor inverts selectivity). That's Evans–Polanyi: the substrate C–H BDE
           difference is abstractor-independent. <strong>So the trunk's advantage is MECHANISTIC</strong>
           (a stable CH₃X decouples activation from over-functionalization; soft sulfur thermodynamics),
           not a "magic catalyst". Product over-functionalization: Cl worse than Br (the bromine route is
           cleaner). <em>(def2-SVP sign-confirmation running on AWS.)</em></p>
           <p><strong>Honest caveats:</strong> thermodynamics = <em>capability</em>, not yield
           (kinetics/coking/separation decide the rest); rigid-scan absolutes are crude (signs are
           robust). The bottleneck shifts to the <em>mediator loop</em> (Cl/HX recycle, H₂S management)
           and CH₃X → C₂ coupling selectivity — the next real AWS calculation.</p>`
    },
    table: {
      title: { ru: "Карта маршрутов: макс. равновесная конверсия CH₄ и класс стены",
               en: "Route map: max equilibrium CH₄ conversion and wall class" },
      head: { ru: ["маршрут", "макс. конв.", "класс"], en: ["route", "max conv.", "wall class"] },
      rows: [
        ["★ Cl₂ → C₂H₆ (галоген/halogen)", "~100% @ 400 °C", "no-COx (мех./mech.)"],
        ["★ Br₂ → C₂H₆ (галоген/halogen)", "~93% @ 400 °C", "no-COx (чище/cleaner)"],
        ["★ S₂ → C₂H₄ (сера/sulfur)", "~61% @ 1400 °C", "no-COx (мягкий/soft)"],
        ["O₂ → C₂ (OCM)", "~100% термо / thermo", "COx wall → ~28% выход/yield"],
        ["non-ox / MDA", "~51–61% @ 1400 °C", "equilibrium-limited"]
      ]
    },
    figures: [
      { src: "assets/route_landscape.png",
        caption: { ru: "Ландшафт маршрутов: зелёные (сера/галоген) обходят COx-стену и термодинамически образуют магистраль; галоген доминирует уже при низкой T.",
                   en: "Route landscape: green (sulfur/halogen) escape the COx wall and form the thermodynamic trunk; halogen dominates already at low T." } }
    ]
  },

  /* ---------------------------------------------------------------- */
  {
    id: "cycles",
    date: "2026-06-06",
    stage: { ru: "Этап 7", en: "Stage 7" },
    title: { ru: "Замыкается ли цикл? Бром выходит вперёд",
             en: "Does the cycle close? Bromine pulls ahead" },
    simple: {
      ru: `<p>Чтобы медиаторный путь был <strong>настоящим процессом</strong>, медиатор
           (хлор/бром/серу) нужно <em>возвращать в цикл</em> — иначе ты его просто расходуешь.
           Проверили термодинамику полного цикла:</p>
           <ul>
             <li><strong>Хлор:</strong> сшивка метана отличная, но вернуть Cl₂ (реакция Дикона)
                 выше ~600 °C <em>застревает</em> — это известная беда хлора.</li>
             <li><strong>Бром:</strong> <em>оба</em> шага идут «под горку» — <strong>цикл
                 замыкается сам.</strong></li>
             <li><strong>Сера:</strong> регенерация тривиальна, но сама сшивка слабая (равновесие).</li>
           </ul>
           <p>Итог: самая чистая самозамыкающаяся магистраль — <strong>БРОМ</strong>. Он соединяет
           метан, не жжёт продукт (нет CO₂), сам себя регенерирует, а суммарная реакция —
           <em>такая же, как у O₂-OCM, но БЕЗ стены COx</em>. Плюс бром меньше «пережёвывает» свой
           продукт, чем хлор (Этап 6). Значит, <strong>бром-медиаторное окислительное сочетание —
           расчётная магистраль</strong>.</p>`,
      en: `<p>For a mediated route to be a <strong>real process</strong>, the mediator
           (chlorine/bromine/sulfur) must be <em>recycled</em> — otherwise you just consume it. We
           checked the full-cycle thermodynamics:</p>
           <ul>
             <li><strong>Chlorine:</strong> methane coupling is excellent, but regenerating Cl₂ (the
                 Deacon step) <em>stalls</em> above ~600 °C — the well-known chlorine problem.</li>
             <li><strong>Bromine:</strong> <em>both</em> steps run downhill — <strong>the cycle
                 closes by itself.</strong></li>
             <li><strong>Sulfur:</strong> regeneration is trivial, but the coupling itself is weak
                 (equilibrium-limited).</li>
           </ul>
           <p>Bottom line: the cleanest self-closing trunk is <strong>BROMINE</strong>. It couples
           methane, doesn't burn the product (no CO₂), regenerates itself, and the net reaction is
           <em>the same as O₂-OCM but WITHOUT the COx wall</em>. Bromine also over-reacts with its
           own product less than chlorine (Stage 6). So <strong>bromine-mediated oxidative coupling is
           the computable trunk</strong>.</p>`
    },
    tech: {
      ru: `<p>ΔG(T) шагов цикла (NIST термохимия, конст. ΔH/ΔS), при 700 °C:</p>
           <p><strong>Cl:</strong> сшивка 2CH₄+Cl₂→C₂H₆+2HCl ΔG=−127; Дикон 2HCl+½O₂→Cl₂+H₂O
           ΔG=+5.6 (пересекает 0 около 600 °C → <em>бутылочное горло</em>). <strong>Br:</strong>
           сшивка −47; Дикон 2HBr+½O₂→Br₂+H₂O ΔG=−74 → <strong>цикл замыкается</strong>.
           <strong>S:</strong> сшивка +2.5 (равновесие); Клаус 2H₂S+O₂→S₂+2H₂O ΔG=−303 → замыкается,
           но лимитирует сшивка. Суммарно для всех: <code>2CH₄+½O₂→C₂H₆+H₂O</code> (или →C₂H₄+H₂O для
           S) — <strong>та же net-реакция, что у OCM, но переокисление обходится медиатором</strong>.</p>
           <p><strong>Победитель — бром:</strong> единственный, у кого ОБА шага термодинамически
           благоприятны (цикл self-regenerating) И минимальное переокисление продукта (Этап 6:
           Br −30 vs Cl −34). Это согласуется с реальными бром-медиаторными процессами конверсии
           метана.</p>
           <p><strong>Честные оговорки:</strong> только термодинамика (конст. ΔH/ΔS); реальную
           селективность и скорость решают кинетика сшивки CH₃•/CH₃Br → C₂ против ре-бромирования
           (концентрационно-зависимо), обращение с HBr/H₂S и коррозия. Это уже инженерия процесса.</p>`,
      en: `<p>ΔG(T) of cycle steps (NIST thermochem, constant ΔH/ΔS), at 700 °C:</p>
           <p><strong>Cl:</strong> coupling 2CH₄+Cl₂→C₂H₆+2HCl ΔG=−127; Deacon 2HCl+½O₂→Cl₂+H₂O
           ΔG=+5.6 (crosses 0 near 600 °C → <em>bottleneck</em>). <strong>Br:</strong> coupling −47;
           Deacon 2HBr+½O₂→Br₂+H₂O ΔG=−74 → <strong>cycle closes</strong>. <strong>S:</strong>
           coupling +2.5 (equilibrium); Claus 2H₂S+O₂→S₂+2H₂O ΔG=−303 → closes, but coupling is the
           limit. Net for all: <code>2CH₄+½O₂→C₂H₆+H₂O</code> (or →C₂H₄+H₂O for S) — <strong>the same
           net reaction as OCM, but over-oxidation is bypassed by the mediator</strong>.</p>
           <p><strong>Winner — bromine:</strong> the only one where BOTH steps are thermodynamically
           favourable (self-regenerating cycle) AND product over-functionalization is lowest (Stage 6:
           Br −30 vs Cl −34). This matches real bromine-mediated methane-conversion routes.</p>
           <p><strong>Honest caveats:</strong> thermodynamics only (constant ΔH/ΔS); real selectivity
           and rate are set by the kinetics of CH₃•/CH₃Br → C₂ coupling vs re-bromination
           (concentration-dependent), HBr/H₂S handling, and corrosion. That's process engineering.</p>`
    },
    table: {
      title: { ru: "Термодинамика цикла (ΔG, кДж/моль, 700 °C)", en: "Cycle thermodynamics (ΔG, kJ/mol, 700 °C)" },
      head: { ru: ["цикл", "сшивка", "регенерация", "вердикт"], en: ["cycle", "coupling", "recycle", "verdict"] },
      rows: [
        ["★ Br (бром/bromine)", "−47", "−74 (Deacon)", "ЦИКЛ ЗАМЫКАЕТСЯ / CLOSES ✅"],
        ["Cl (хлор/chlorine)", "−127", "+5.6 (Deacon)", "регенерация-горло / recycle bottleneck"],
        ["S (сера/sulfur)", "+2.5", "−303 (Claus)", "сшивка-горло / coupling-limited"]
      ]
    },
    figures: [
      { src: "assets/cycle_recycle.png",
        caption: { ru: "ΔG регенерации медиатора vs T: бром (оранж.) под горку везде → цикл замыкается; хлор (синий) выше ~600°C застревает; сера (Клаус) тривиально замыкается.",
                   en: "ΔG of mediator recycle vs T: bromine (orange) downhill everywhere → cycle closes; chlorine (blue) stalls above ~600 °C; sulfur (Claus) closes trivially." } }
    ]
  },

  /* ---------------------------------------------------------------- */
  {
    id: "ethylene",
    date: "2026-06-06",
    stage: { ru: "Этап 8", en: "Stage 8" },
    title: { ru: "Этилен — слабое звено: медиаторы едят двойную связь",
             en: "Ethylene is the weak link: mediators eat the double bond" },
    simple: {
      ru: `<p>Нашли, <strong>почему «один шаг до этилена» так труден.</strong> Слабое место
           этилена — не его C–H связи (они как раз <em>прочные</em>), а его <strong>двойная
           связь</strong>. И хлор/бром/сера радостно <em>атакуют именно её</em> (присоединяются),
           разрушая этилен.</p>
           <p>Отсюда честный вывод про бром: он чисто делает <strong>ЭТАН</strong> (у этана нет
           двойной связи → он в безопасности; цикл замыкается, COx-стены нет). Но как только
           появляется <strong>этилен — бром его съедает.</strong> Значит, бромная магистраль =
           метан→этан (отлично, обходит стену), а этилен требует <strong>отдельной
           бром-свободной стадии</strong> (дегидрирования этана).</p>
           <p><strong>«Один шаг до этилена» блокируется дважды:</strong> в O₂-путях — сгоранием в
           CO₂; в медиаторных — присоединением к двойной связи. <em>Просвет:</em> кислородный
           радикал (OH) двойную связь НЕ предпочитает (+1.7, «щадит») — значит, может существовать
           особый медиатор, который активирует метан, но щадит этилен. Это следующая охота.</p>`,
      en: `<p>We found <strong>why "one step to ethylene" is so hard.</strong> Ethylene's weak spot
           is not its C–H bonds (those are actually <em>strong</em>) but its <strong>double
           bond</strong>. And chlorine/bromine/sulfur happily <em>attack exactly that</em> (add to
           it), destroying the ethylene.</p>
           <p>Honest consequence for bromine: it cleanly makes <strong>ETHANE</strong> (ethane has no
           double bond → it's safe; the cycle closes, no COx wall). But the moment <strong>ethylene
           appears, bromine eats it.</strong> So the bromine trunk = methane→ethane (great, escapes
           the wall), while ethylene needs a <strong>separate bromine-free step</strong> (dehydrogenation
           of ethane).</p>
           <p><strong>"One step to ethylene" is blocked twice:</strong> in O₂ routes by burning to CO₂;
           in mediated routes by addition to the double bond. <em>A hopeful sign:</em> the oxygen
           radical (OH) does NOT prefer the double bond (+1.7, "spares" it) — so a special mediator
           that activates methane yet spares ethylene might exist. That's the next hunt.</p>`
    },
    tech: {
      ru: `<p>Энергии реакций (def2-SVP, оптимизированные минимумы; Эванс–Поляни: ниже ΔH →
           ниже барьер → быстрее). Дескриптор = ΔΔH относительно активации метана.</p>
           <p><strong>Винильный HAT</strong> X•+C₂H₄→XH+C₂H₃•: ΔΔH = <strong>+8.7 для ВСЕХ</strong>
           (X сокращается → это просто разница BDE винил-vs-метил C–H; винильная C–H прочнее → НЕ
           тот канал). <strong>π-присоединение</strong> X•+C₂H₄→X-CH₂CH₂•: ΔΔH = OH <strong>+1.7
           (щадит)</strong>, Cl <strong>−31.3</strong>, Br <strong>−19.5</strong>, SH <strong>−19.6
           (разрушают)</strong>.</p>
           <p><strong>Следствие:</strong> медиаторная сшивка даёт насыщенный продукт (этан) чисто,
           но ненасыщенную цель (этилен) Cl/Br/S съедают π-присоединением → метан→этилен в один шаг
           медиатором заблокирован. Архитектура: бром→этан (цикл, Этап 7) + отдельное бром-свободное
           дегидрирование→этилен. O/OH π-связь щадит, но у O₂-пути своя COx-стена. <em>Оговорка:</em>
           прокси по энергиям реакций (без явного TS); знаки надёжны.</p>`,
      en: `<p>Reaction energies (def2-SVP, optimized minima; Evans–Polanyi: lower ΔH → lower barrier
           → faster). Descriptor = ΔΔH relative to methane activation.</p>
           <p><strong>Vinylic HAT</strong> X•+C₂H₄→XH+C₂H₃•: ΔΔH = <strong>+8.7 for ALL</strong> (X
           cancels → just the vinyl-vs-methyl C–H BDE gap; vinylic C–H is stronger → NOT the channel).
           <strong>π-addition</strong> X•+C₂H₄→X-CH₂CH₂•: ΔΔH = OH <strong>+1.7 (spares)</strong>,
           Cl <strong>−31.3</strong>, Br <strong>−19.5</strong>, SH <strong>−19.6 (destroy)</strong>.</p>
           <p><strong>Implication:</strong> mediated coupling gives the saturated product (ethane)
           cleanly, but the unsaturated target (ethylene) is eaten by Cl/Br/S via π-addition → one-step
           methane→ethylene by a mediator is blocked. Architecture: bromine→ethane (cycle, Stage 7) +
           a separate bromine-free dehydrogenation→ethylene. O/OH spares the π-bond, but the O₂ route
           has its own COx wall. <em>Caveat:</em> reaction-energy proxy (no explicit TS); signs robust.</p>`
    },
    table: {
      title: { ru: "Как медиатор атакует этилен (ΔΔH vs активации метана, ккал/моль)",
               en: "How the mediator attacks ethylene (ΔΔH vs methane activation, kcal/mol)" },
      head: { ru: ["медиатор", "π-присоед. к C=C", "винильный HAT", "этилен"],
              en: ["mediator", "π-add to C=C", "vinylic HAT", "ethylene"] },
      rows: [
        ["OH• (кислород/oxygen)", "+1.7", "+8.7", "ЩАДИТ / SPARED"],
        ["Cl•", "−31.3", "+8.7", "РАЗРУШЕН / DESTROYED"],
        ["Br•", "−19.5", "+8.7", "РАЗРУШЕН / DESTROYED"],
        ["SH• (сера/sulfur)", "−19.6", "+8.7", "РАЗРУШЕН / DESTROYED"]
      ]
    },
    figures: [
      { src: "assets/ethylene_channel.png",
        caption: { ru: "π-присоединение (красное) разрушает этилен у Cl/Br/S; винильный HAT (синее, +9) — не тот канал; OH этилен щадит.",
                   en: "π-addition (red) destroys ethylene for Cl/Br/S; vinylic HAT (blue, +9) is the wrong channel; OH spares ethylene." } }
    ]
  },

  /* ---------------------------------------------------------------- */
  {
    id: "hunt",
    date: "2026-06-06",
    stage: { ru: "Этап 9", en: "Stage 9" },
    title: { ru: "Охота на медиатор: кто активирует метан, но щадит этилен?",
             en: "The mediator hunt: who activates methane yet spares ethylene?" },
    simple: {
      ru: `<p>Святой Грааль одного шага: радикал-«партнёр», который <em>отнимает водород у метана</em>
           (чтобы запустить сочетание), но <strong>НЕ трогает двойную связь этилена</strong> (чтобы не
           разрушить продукт). Прогнали кислородные/азотные/фторные радикалы.</p>
           <p>Галогены, сера, фтор — все <strong>едят</strong> двойную связь. Но <strong>N/O-радикалы
           (NH₂, OH, CF₃) её ЩАДЯТ!</strong> Звезда — <strong>NH₂</strong> (азотный радикал): активирует
           метан <em>и</em> щадит этилен, и (в отличие от кислорода) у него нет пути «сгорания в CO₂».</p>
           <p>→ <strong>Азот-медиаторный путь — свежая вычисляемая гипотеза для одного шага.</strong>
           Честные оговорки: NH₂ активирует метан лишь слабо (медленно); азот-радикальные медиаторы и их
           рецикл — экзотика, неизучено; кислород (OH) щадит π-связь, но жжёт продукт <em>другими</em>
           путями (та самая COx-стена). Это <strong>лид, а не решение</strong> — новое направление,
           которое движок нащупал сам.</p>`,
      en: `<p>The holy grail of the one-step route: a "partner" radical that <em>takes a hydrogen from
           methane</em> (to start coupling) but <strong>does NOT touch ethylene's double bond</strong>
           (so it won't destroy the product). We screened oxygen/nitrogen/fluorine radicals.</p>
           <p>Halogens, sulfur, fluorine all <strong>eat</strong> the double bond. But <strong>N/O
           radicals (NH₂, OH, CF₃) SPARE it!</strong> The standout is <strong>NH₂</strong> (a nitrogen
           radical): it activates methane <em>and</em> spares ethylene, and (unlike oxygen) has no
           "burn-to-CO₂" path.</p>
           <p>→ <strong>A nitrogen-mediated route is a fresh, computable hypothesis for the one-step
           goal.</strong> Honest caveats: NH₂ activates methane only weakly (slow); nitrogen-radical
           mediators and their recycle are exotic/unexplored; oxygen (OH) spares the π-bond but burns the
           product by <em>other</em> paths (the COx wall). This is a <strong>lead, not a solution</strong>
           — a new direction the engine surfaced on its own.</p>`
    },
    tech: {
      ru: `<p>Энергии реакций (def2-SVP, опт. минимумы). Зона-энэйблер (активация &lt;~3 И ddH_π&gt;0):
           <strong>OH</strong> (−8.4, +1.7), <strong>NH₂</strong> (−0.8, +1.4), <strong>CF₃</strong>
           (+2.1, +2.6). Едят этилен (ddH_π&lt;0): Cl (−31.3), Br (−19.5), SH (−19.6), F (−1.3), OOH (−0.2).
           CH₃O — на грани (0.0).</p>
           <p><strong>Лид — NH₂:</strong> единственный с благоприятной активацией метана, щадящий π-связь,
           и БЕЗ кислорода → нет COx-пути глубокого окисления. OH тоже щадит π, но O жжёт C₂ другими путями
           (COx-стена — не ловится одним π-дескриптором). CF₃ щадит, но активация вверх (медленно). F
           слишком агрессивен.</p>
           <p><strong>Гипотеза (свежая, вычисляемая):</strong> азот-центрированная радикальная медиация
           для метан→этилен в один шаг. Нужно проверить TS-уровнем (явные барьеры), продумать источник/
           рецикл NH₂-медиатора, и убедиться в отсутствии N-специфичного переокисления. <em>Оговорка:</em>
           прокси по энергиям реакций; знаки/тренды надёжны, абсолюты — нет.</p>`,
      en: `<p>Reaction energies (def2-SVP, optimized minima). Enabler zone (activation &lt;~3 AND ddH_π&gt;0):
           <strong>OH</strong> (−8.4, +1.7), <strong>NH₂</strong> (−0.8, +1.4), <strong>CF₃</strong>
           (+2.1, +2.6). Eat ethylene (ddH_π&lt;0): Cl (−31.3), Br (−19.5), SH (−19.6), F (−1.3), OOH (−0.2).
           CH₃O is borderline (0.0).</p>
           <p><strong>The lead — NH₂:</strong> the only one with favourable methane activation that spares
           the π-bond AND has no oxygen → no COx deep-oxidation path. OH also spares π, but O burns C₂ by
           other paths (the COx wall — not captured by a single π descriptor). CF₃ spares but activation is
           uphill (slow). F is too aggressive.</p>
           <p><strong>Hypothesis (fresh, computable):</strong> nitrogen-centred radical mediation for
           methane→ethylene in one step. Needs TS-level verification (explicit barriers), a source/recycle
           for the NH₂ mediator, and a check for N-specific over-oxidation. <em>Caveat:</em> reaction-energy
           proxy; signs/trends robust, absolutes not.</p>`
    },
    table: {
      title: { ru: "Карта медиаторов (ккал/моль)", en: "Mediator map (kcal/mol)" },
      head: { ru: ["медиатор", "активация CH₄", "щадит этилен (ddH_π)", "зона"],
              en: ["mediator", "CH₄ activation", "spares ethylene (ddH_π)", "zone"] },
      rows: [
        ["★ NH₂• (азот/nitrogen)", "−0.8", "+1.4", "ЭНЭЙБЛЕР / ENABLER (no COx)"],
        ["OH• (кислород/oxygen)", "−8.4", "+1.7", "щадит π, но COx-стена / spares π but COx wall"],
        ["CF₃•", "+2.1", "+2.6", "щадит, активация вверх / spares, slow"],
        ["Cl / Br / SH", "плохо/poor", "−20…−31", "ЕДЯТ / EAT C=C"],
        ["F•", "−21.3", "−1.3", "агрессивен, ест / aggressive, eats"]
      ]
    },
    figures: [
      { src: "assets/mediator_map.png",
        caption: { ru: "Карта медиаторов: верх-лево = активирует CH₄ + щадит C₂H₄ (зона одного шага). NH₂/OH/CF₃ щадят, галогены/S/F едят. NH₂ — настоящий лид (без COx).",
                   en: "Mediator map: upper-left = activates CH₄ + spares C₂H₄ (one-step zone). NH₂/OH/CF₃ spare, halogens/S/F eat. NH₂ is the genuine lead (no COx)." } }
    ]
  },

  /* ---------------------------------------------------------------- */
  {
    id: "verdict",
    date: "2026-06-06",
    stage: { ru: "Этап 10", en: "Stage 10" },
    title: { ru: "Честный вердикт: карта всех путей и где магистраль",
             en: "The honest verdict: the map of every path, and where the trunk is" },
    simple: {
      ru: `<p>Весь путь на одной картинке. Мечта была: метан → этилен в ОДИН шаг при 70%. Честный итог:</p>
           <ul>
             <li><strong>O₂-OCM (28% стена), электроника, сито, конфайнмент — стену НЕ ломают</strong>
                 (артефакты мы ловили контролями).</li>
             <li><strong>Революционный выход из C₂-стены — медиаторный цикл:</strong> бром соединяет метан
                 в <em>этан</em>, обходит «сгорание» и сам регенерируется (магистраль).</li>
             <li><strong>Но один шаг именно до ЭТИЛЕНА всё ещё заблокирован</strong> (O₂ — сгорание;
                 медиаторы — атака двойной связи). Этилен = бром→этан + отдельное дегидрирование (два шага).</li>
             <li><strong>Свежий лид для одного шага:</strong> азот-медиатор (NH₂-тип), который щадит двойную
                 связь — новое, вычисляемое, но непроверенное направление.</li>
           </ul>
           <p><strong>Итог честно:</strong> вычисляемая революционная магистраль <em>есть</em> (бром→этан,
           обходит 40-летнюю стену); «один шаг до этилена при 70%» остаётся фронтиром — но теперь с
           конкретным новым направлением (азот-медиаторы), которое стоит копать.</p>`,
      en: `<p>The whole journey in one picture. The dream was methane → ethylene in ONE step at 70%. The
           honest verdict:</p>
           <ul>
             <li><strong>O₂-OCM (28% wall), electronics, sieving, confinement — do NOT break the wall</strong>
                 (we caught the artefacts with controls).</li>
             <li><strong>The revolutionary escape from the C₂ wall is a mediated cycle:</strong> bromine
                 couples methane to <em>ethane</em>, routes around the burning wall, and self-regenerates
                 (the trunk).</li>
             <li><strong>But one step to ETHYLENE specifically is still blocked</strong> (O₂ — burning;
                 mediators — double-bond attack). Ethylene = bromine→ethane + a separate dehydrogenation
                 (two steps).</li>
             <li><strong>A fresh lead for the one-step goal:</strong> a nitrogen mediator (NH₂-type) that
                 spares the double bond — a new, computable, but unproven direction.</li>
           </ul>
           <p><strong>Honest bottom line:</strong> a computable revolutionary trunk <em>exists</em>
           (bromine→ethane, escapes the 40-year wall); "one step to ethylene at 70%" remains the frontier —
           but now with a concrete new direction (nitrogen mediators) worth chasing.</p>`
    },
    tech: {
      ru: `<p>Сводка по всем путям (лучший честно-достижимый выход за проход / стена, что его держит):
           электроника ~29% (BEP-стена); сито ~28% (CH₄≈C₂H₄ по размеру); конфайнмент — артефакт (не
           подтверждён контролем на def2-SVP); неокислительный ~50% (равновесие); сера ~55%
           (сшивка-равновесие); <strong>бром-медиатор → этан: цикл замыкается (Дикон −74), COx-стены нет
           (магистраль)</strong>; тот же путь → этилен — блок π-присоединением; <strong>NH₂-медиатор —
           лид для одного шага</strong> (щадит π, без COx; не доказан).</p>
           <p>Все числа — из реальных расчётов в репозитории (термодинамика NIST + газофазная квант-химия
           CASCI/UKS), с гейтами честности и контролями (дальний воротник, парность, сходимость; артефакты
           помечены/отброшены). Квантового превосходства не заявляется. Запуски на AWS, суммарно ~$3.</p>
           <p><strong>Следующий настоящий расчёт:</strong> TS-уровень для NH₂-медиации (явные барьеры
           активации/π-щадения) + источник и рецикл азот-медиатора + проверка N-переокисления; параллельно
           — кинетика сшивки vs ре-бромирования для бром→этан (концентрационная селективность).</p>`,
      en: `<p>Summary across all paths (best honestly-achievable single-pass yield / the wall that stops
           it): electronics ~29% (BEP wall); sieving ~28% (CH₄≈C₂H₄ in size); confinement — artefact (not
           confirmed by the def2-SVP control); non-oxidative ~50% (equilibrium); sulfur ~55%
           (coupling-equilibrium); <strong>bromine mediator → ethane: cycle closes (Deacon −74), no COx wall
           (the trunk)</strong>; the same route → ethylene is blocked by π-addition; <strong>NH₂ mediator —
           the lead for the one-step goal</strong> (spares π, no COx; unproven).</p>
           <p>Every number comes from a real computation in the repo (NIST thermodynamics + gas-phase
           quantum chemistry CASCI/UKS), with honesty gates and controls (far-collar, parity, convergence;
           artefacts flagged/discarded). No quantum advantage is claimed. Runs on AWS, ~$3 total.</p>
           <p><strong>The next real calculation:</strong> TS-level study of NH₂ mediation (explicit
           activation / π-sparing barriers) + a source and recycle for the nitrogen mediator + a check for
           N-over-oxidation; in parallel, the kinetics of coupling vs re-bromination for bromine→ethane
           (concentration-driven selectivity).</p>`
    },
    figures: [
      { src: "assets/synthesis_map.png",
        caption: { ru: "Честная карта всех путей: стену ломает медиаторный цикл (бром→этан); один шаг до этилена пока заблокирован; азот-медиатор (NH₂) — свежий лид.",
                   en: "The honest map of every path: the wall is broken by a mediated cycle (bromine→ethane); one-step ethylene is still blocked; the nitrogen mediator (NH₂) is the fresh lead." } }
    ]
  },

  /* ---------------------------------------------------------------- */
  {
    id: "nfamily",
    date: "2026-06-06",
    stage: { ru: "Этап 11", en: "Stage 11" },
    title: { ru: "Лид под давлением: у азота есть компромисс",
             en: "The lead under pressure: nitrogen has a trade-off" },
    simple: {
      ru: `<p>Стресс-тест азотного лида. Вопрос: «азот щадит этилен» — это общее свойство или
           случайность именно NH₂? Проверили электрон-богатый (NMe₂) и электрон-бедный (NF₂)
           азотные радикалы.</p>
           <p><strong>Хорошая новость:</strong> щадение этилена ОБОБЩАЕТСЯ — и NH₂, и NMe₂ берегут
           двойную связь. Это реальное свойство азота → дизайн-правило.</p>
           <p><strong>Засада — компромисс:</strong> делаешь азот «мягче» (богаче электронами, NMe₂)
           → ещё лучше щадит этилен, но становится <em>слишком слабым, чтобы оторвать водород у
           метана</em> (не активирует). Делаешь «жёстче» (беднее, NF₂) → вообще метан не берёт.
           <strong>NH₂ — единственная «сладкая точка»</strong>: ровно настолько реактивен, чтобы
           взять метан и при этом пощадить этилен — и даже он на грани (медленно).</p>
           <p>Честно: лид настоящий, но живёт в <strong>узком окне</strong>, и лучший кандидат
           кинетически слабоват. Следующее: откуда брать NH₂ и можно ли его вернуть в цикл
           (жизнеспособность), и проверка барьеров на уровне переходных состояний.</p>`,
      en: `<p>Stress-testing the nitrogen lead. Question: is "nitrogen spares ethylene" a general
           property, or special to NH₂? We tested an electron-rich (NMe₂) and an electron-poor (NF₂)
           nitrogen radical.</p>
           <p><strong>Good news:</strong> the π-sparing GENERALIZES — both NH₂ and NMe₂ spare the
           double bond. It's a real nitrogen property → a design rule.</p>
           <p><strong>The catch — a trade-off:</strong> make the nitrogen "softer" (electron-rich,
           NMe₂) → it spares ethylene even more, but becomes <em>too weak to take methane's
           hydrogen</em> (no activation). Make it "harder" (electron-poor, NF₂) → it can't activate
           methane at all. <strong>NH₂ is the only "sweet spot"</strong>: just reactive enough to grab
           methane while sparing ethylene — and even it is borderline (slow).</p>
           <p>Honest: the lead is real but lives in a <strong>narrow window</strong>, and the best
           candidate is kinetically weak. Next: where NH₂ comes from and whether it can be recycled
           (viability), plus a transition-state-level barrier check.</p>`
    },
    tech: {
      ru: `<p>Активация / ddH_π (def2-SVP, опт. минимумы): <strong>NH₂</strong> (−0.8, +1.4 щадит),
           <strong>NMe₂</strong> (+12.3, +2.2 щадит), <strong>NF₂</strong> (+67.3, −67.8 — но это
           доминируется провальной активацией, не «поеданием» этилена; π-add у NF₂ сам по себе ~0).</p>
           <p><strong>Тренд:</strong> π-щадение общее для амино-радикалов; активация и π-реактивность
           обе растут с H-сродством азотного радикала → <strong>антикорреляция «активация ↔ щадение»</strong>.
           Богатый электронами N (NMe₂) — слабый отрыватель H (плохая активация) но отлично щадит;
           бедный (NF₂) — N–H связь в продукте слаба → активация невозможна. NH₂ — единственное
           пересечение «активирует И щадит», и то у самой границы (−0.8 → кинетически медленно).</p>
           <p><strong>Оговорки:</strong> прокси по энергиям реакций (без явного TS — барьеры могут
           сдвинуть окно); нужна термодинамика источника/рецикла NH₂ (N–H BDE в NH₃ ~107 ккал/моль —
           источник дорог) и TS-уровень. Лид реален как направление, но с жёсткими ограничениями.</p>`,
      en: `<p>Activation / ddH_π (def2-SVP, optimized minima): <strong>NH₂</strong> (−0.8, +1.4 spares),
           <strong>NMe₂</strong> (+12.3, +2.2 spares), <strong>NF₂</strong> (+67.3, −67.8 — but this is
           dominated by the failed activation, not by "eating" ethylene; NF₂'s π-add alone is ~0).</p>
           <p><strong>Trend:</strong> π-sparing is general for amine radicals; activation and
           π-reactivity both scale with the N-radical's H-affinity → an <strong>activation↔sparing
           anticorrelation</strong>. An electron-rich N (NMe₂) is a weak H-abstractor (poor activation)
           but spares well; an electron-poor N (NF₂) has a weak product N–H → activation impossible. NH₂
           is the only intersection of "activates AND spares", and it's at the very edge (−0.8 →
           kinetically slow).</p>
           <p><strong>Caveats:</strong> reaction-energy proxy (no explicit TS — barriers could shift the
           window); needs the thermodynamics of an NH₂ source/recycle (the N–H BDE in NH₃ ~107 kcal/mol
           — an expensive source) and a TS-level check. The lead is real as a direction, but tightly
           constrained.</p>`
    },
    table: {
      title: { ru: "N-семейство: активация vs щадение (ккал/моль)", en: "N-family: activation vs sparing (kcal/mol)" },
      head: { ru: ["N-радикал", "активация CH₄", "ddH_π (щадит?)", "итог"], en: ["N-radical", "CH₄ activation", "ddH_π (spares?)", "verdict"] },
      rows: [
        ["★ NH₂", "−0.8", "+1.4 ✓", "сладкая точка (на грани) / sweet spot (borderline)"],
        ["NMe₂ (e-rich)", "+12.3", "+2.2 ✓", "щадит, но не активирует / spares, no activation"],
        ["NF₂ (e-poor)", "+67.3", "n/a", "не активирует метан / can't activate methane"]
      ]
    },
    figures: [
      { src: "assets/nfamily_tradeoff.png",
        caption: { ru: "Компромисс азота: все амино-радикалы щадят этилен, но активация падает при настройке на щадение. Только NH₂ попадает в окно «активирует И щадит».",
                   en: "Nitrogen's trade-off: all amine radicals spare ethylene, but activation collapses as you tune for sparing. Only NH₂ lands in the 'activates AND spares' window." } }
    ]
  },

  /* ---------------------------------------------------------------- */
  {
    id: "viability",
    date: "2026-06-06",
    stage: { ru: "Этап 12", en: "Stage 12" },
    title: { ru: "Можно ли вернуть азот в цикл? Нет — лид не замыкается",
             en: "Can nitrogen be recycled? No — the lead doesn't close" },
    simple: {
      ru: `<p>Финальная проверка азотного лида. Как брому нужен шаг Дикона, чтобы вернуть Br₂,
           так медиатор NH₂ должен регенерироваться из NH₃ — иначе ты его просто расходуешь.
           Посчитали:</p>
           <ul>
             <li>вернуть NH₂ из NH₃ — <strong>сильно «в гору»</strong> (+59 гомолизом, +64 даже с
                 кислородом): связь N–H слишком прочная;</li>
             <li>единственный выгодный путь (через OH-радикал, −8) <strong>требует кислорода</strong>
                 — а тот возвращает «сгорание» (COx-проблему).</li>
           </ul>
           <p><strong>Итог:</strong> азотный лид <em>щадит</em> этилен (здорово), но <strong>не
           замыкается в цикл</strong> — NH₂ нельзя дёшево регенерировать без кислорода. В отличие
           от брома, чей цикл <em>замыкается</em> (Дикон −18). Честное закрытие: самый красивый
           лид на «один шаг», доведённый до проверки жизнеспособности, не выживает.</p>
           <p>Так что бром→этан остаётся единственной <em>жизнеспособной</em> магистралью; «один
           шаг до этилена 70%» честно остаётся непокорённым по всем просчитанным путям. Это и есть
           полная, честная карта — с одной реальной дверью (бром→этан) и ясно очерченными стенами.</p>`,
      en: `<p>The final test of the nitrogen lead. Just as bromine needs its Deacon step to
           regenerate Br₂, the NH₂ mediator must be regenerated from NH₃ — else you just consume it.
           Computed:</p>
           <ul>
             <li>regenerating NH₂ from NH₃ is strongly <strong>uphill</strong> (+59 by homolysis,
                 +64 even with oxygen): the N–H bond is too strong;</li>
             <li>the only favourable route (via an OH radical, −8) <strong>needs oxygen</strong> —
                 which brings back the "burning" (COx) problem.</li>
           </ul>
           <p><strong>Verdict:</strong> the nitrogen lead <em>spares</em> ethylene (great) but
           <strong>doesn't close into a cycle</strong> — NH₂ can't be cheaply regenerated without
           oxygen. Unlike bromine, whose cycle <em>closes</em> (Deacon −18). An honest closure: the
           most beautiful one-step lead, taken all the way to a viability check, does not survive.</p>
           <p>So bromine→ethane remains the one <em>viable</em> trunk; "one step to ethylene at 70%"
           honestly stays unbeaten across every path computed. That is the complete, honest map — one
           real door (bromine→ethane) and clearly-drawn walls.</p>`
    },
    tech: {
      ru: `<p>Энергии регенерации NH₂ (def2-SVP/PBE0, ΔE): NH₃→NH₂•+½H₂ <strong>+59.5</strong>;
           2NH₃+½O₂→2NH₂•+H₂O <strong>+64.4</strong>; NH₃+OH•→NH₂•+H₂O <strong>−7.6</strong>.
           Прочная N–H связь (BDE ~107 ккал/моль) делает регенерацию эндотермичной; O₂-ассист не
           спасает; только OH-ассист выгоден, но возвращает O-радикальную химию (COx). Сравнение:
           бромный Дикон (Этап 7) замыкается (ΔG −74 кДж/моль ≈ −18 ккал/моль).</p>
           <p><strong>Вывод:</strong> N-медиатор = щадит π-связь, но не рециклируется → не
           жизнеспособный цикл. Полная карта: единственная замыкающаяся медиаторная магистраль —
           бром→этан; один шаг до этилена заблокирован везде (O₂ — сгорание; Cl/Br/S — π-присоединение;
           N — не рециклируется). <em>Оговорка:</em> термодинамика; не оценены неравновесные источники
           N-радикалов (плазма/фото) — спекулятивный обход, не считался.</p>`,
      en: `<p>NH₂ regeneration energies (def2-SVP/PBE0, ΔE): NH₃→NH₂•+½H₂ <strong>+59.5</strong>;
           2NH₃+½O₂→2NH₂•+H₂O <strong>+64.4</strong>; NH₃+OH•→NH₂•+H₂O <strong>−7.6</strong>. The
           strong N–H bond (BDE ~107 kcal/mol) makes regeneration endothermic; O₂-assist doesn't
           rescue it; only OH-assist is favourable but reintroduces O-radical chemistry (COx).
           Contrast: the bromine Deacon step (Stage 7) closes (ΔG −74 kJ/mol ≈ −18 kcal/mol).</p>
           <p><strong>Conclusion:</strong> the N-mediator spares the π-bond but doesn't recycle → not a
           viable cycle. The full map: the only self-closing mediated trunk is bromine→ethane; one-step
           ethylene is blocked everywhere (O₂ — burning; Cl/Br/S — π-addition; N — no recycle).
           <em>Caveat:</em> thermodynamics; non-equilibrium N-radical sources (plasma/photo) not
           assessed — a speculative bypass, not computed.</p>`
    },
    table: {
      title: { ru: "Регенерация NH₂-медиатора (ΔE, ккал/моль)", en: "NH₂-mediator regeneration (ΔE, kcal/mol)" },
      head: { ru: ["путь регенерации", "ΔE", "вердикт"], en: ["regeneration route", "ΔE", "verdict"] },
      rows: [
        ["NH₃ → NH₂• + ½H₂", "+59.5", "в гору — не идёт / uphill — fails"],
        ["2NH₃ + ½O₂ → 2NH₂• + H₂O", "+64.4", "в гору даже с O₂ / fails even with O₂"],
        ["NH₃ + OH• → NH₂• + H₂O", "−7.6", "ок, но нужен O→COx / ok, but needs O→COx"],
        ["(сравн.) Br Дикон / Br Deacon", "≈ −18", "ЗАМЫКАЕТСЯ / CLOSES ✅"]
      ]
    },
    figures: []
  },

  /* ---------------------------------------------------------------- */
  {
    id: "real-center",
    date: "2026-06-08",
    accent: "#f6ad55",
    stage: { ru: "Этап 13", en: "Stage 13" },
    title: { ru: "От учебного Li/MgO к настоящему Mn/Na/W: где DFT обязан врать",
             en: "From the training Li/MgO to a real Mn/Na/W centre: where DFT must lie" },
    simple: {
      ru: `<p>До сих пор реальные расчёты мы вели на «учебном» катализаторе <strong>Li/MgO</strong> (классика Лансфорда) — на нём удобно отлаживать движок. Но промышленный OCM-катализатор другой: <strong>Mn–Na₂WO₄/SiO₂</strong>. Главное обещание OCM — честно показать, где обычные методы (DFT) врут, — надо проверять именно на нём.</p>
           <p>Мы собрали минимальную модель его «рабочего кислорода» — фрагмент <strong>Mn=O</strong>, который отрывает водород у метана — и задали прямой вопрос: <strong>можно ли тут доверять обычному DFT?</strong></p>
           <p>Ответ виден прямо в «числах заполнения» орбиталей. У послушной молекулы они равны 0 или 2 («пусто» / «занято»). Здесь же мы <strong>посчитали</strong> и получили дробные: <strong>1.70, 1.27, 0.73, 0.30</strong> — электроны «размазаны» по нескольким состояниям сразу. <strong>Один-единственный «снимок» электронов, на котором стоит обычный DFT, такое описать не может</strong> — он обязан ошибаться. Это и есть наша фишка отличия от IBM/Google: мы показываем не лозунгом, а <em>измеренным числом</em>, где и почему классика врёт.</p>
           <p><strong>Честный итог по «стене 28%»:</strong> мы доказали, что задача мультиреференсна и DFT здесь ненадёжен — это качественный error-bar. А точную величину поправки «на коленке» посчитать нельзя: упрощённые прикидки дали явные артефакты (мы их показываем, а не прячем), поэтому аккуратное число — это следующий шаг (более тяжёлый расчёт / квантовый компьютер).</p>
           <p><strong>И мы отметили, где классике конец:</strong> наша маленькая модель компьютеру пока «по зубам» — никакого квантового превосходства сегодня, говорим прямо. Но настоящий центр требует гораздо большего расчёта, и точная классика упирается в стену около <strong>40 кубитов</strong> — вот с этого порога реально нужен квантовый компьютер. Фронтирный VQE мы гоняем в облаке AWS.</p>`,
      en: `<p>Until now our real calculations ran on the "training" catalyst <strong>Li/MgO</strong> (the classic Lunsford site) — handy for debugging the engine. But the industrial OCM catalyst is different: <strong>Mn–Na₂WO₄/SiO₂</strong>. OCM's core promise — to show honestly where ordinary methods (DFT) fail — must be tested on <em>it</em>.</p>
           <p>We built a minimal model of its "working oxygen" — an <strong>Mn=O</strong> fragment that abstracts an H from methane — and asked a blunt question: <strong>can ordinary DFT be trusted here?</strong></p>
           <p>The answer is right there in the orbital "occupation numbers". For a well-behaved molecule they are 0 or 2 ("empty"/"filled"). Here we <strong>computed</strong> them and got fractional values: <strong>1.70, 1.27, 0.73, 0.30</strong> — the electrons are "smeared" across several states at once. <strong>A single electronic snapshot, which ordinary DFT relies on, cannot represent this</strong> — it must err. This is our differentiation from IBM/Google: we show, with a <em>measured number</em> rather than a slogan, where and why the classical method lies.</p>
           <p><strong>Honest verdict on the "28% wall":</strong> we proved the problem is multireference and DFT is unreliable here — a qualitative error-bar. The exact size of the correction can't be had on the back of an envelope: our shortcut estimates gave clear artefacts (which we show, not hide), so a trustworthy number is the next step (a heavier calculation / a quantum computer).</p>
           <p><strong>And we marked where classical ends:</strong> our small model is still easy for a computer — no quantum advantage today, said plainly. But the real centre needs a far larger calculation, and exact classical diagonalization hits a wall around <strong>40 qubits</strong> — that is the honest threshold from which a quantum computer is genuinely needed. We run frontier VQE in the AWS cloud.</p>`
    },
    tech: {
      ru: `<p><strong>Центр:</strong> минимальная кластерная модель радикально-кислородного HAT-сайта Mn–Na₂WO₄/SiO₂ — фрагмент <strong>Mn=O</strong> (секстет), def2-SVP, на нашем же валидированном движке (PySCF, ROHF→AVAS→CASCI). <em>НЕ</em> Li/MgO. Знаковая конвенция дескриптора — как в «карте потолка»: <code>ΔΔE‡ = барьер(C₂H₆) − барьер(CH₄)</code>, &lt;0 = доокисляет C₂-продукт быстрее метана (стена BEP), для 70% нужно ≳ +4.</p>
           <p><strong>NOON-диагностика (базис-инвариантная) — ядро результата:</strong> в HAT-точке заполнения натуральных орбиталей AVAS(Mn 3d + O 2p) сильно дробные: <code>[2.00, 1.70, 1.27, 1.00, 1.00, 1.00, 0.73, 0.30]</code>. Тройка ~1.00 — это высокоспиновые однократно занятые орбитали (секстет), <em>не</em> мультиреференс; а вот пары <strong>1.70/0.30 и 1.27/0.73</strong> — настоящая статическая корреляция (≈ один-два «эффективно неспаренных» электрона сверх спина). Однодетерминантный DFT этого воспроизвести не может — измерено, а не заявлено.</p>
           <p><strong>Error-bar — честно про метод (провалы показаны, не сглажены):</strong> качественный вывод твёрдый. Количественную поправку к ΔΔE‡ «в лоб» получить не удалось, и мы показываем почему: (1) жёсткий скан «протягивал» H за равновесие O–H → ложный TS на краю (барьер CH₄ раздут, C₂H₆ ушёл «под горку» — артефакт геометрии); (2) AVAS давал РАЗНЫЕ активные пространства в реагенте и TS (CAS(11e,9o) vs CAS(9e,8o)) → барьер CASCI не определён; (3) фиксированное фронтирное CAS(17e,12o) сходится, но малое CAS теряет динамическую корреляцию — CASCI-барьер CH₄ вышел <strong>94</strong> против <strong>11</strong> ккал/моль у DFT (разрыв ~80 ккал/моль — это <em>недостача динамической корреляции, а не физика селективности</em>). Вывод: малое CASCI — это <strong>диагностика (NOON)</strong>, а не барьер; честная количественная поправка требует сбалансированного многоконфигурационного метода (CASPT2/NEVPT2) или VQE на реальном центре. Это и есть содержание стадии «квантовая поправка».</p>
           <p><strong>Где нужен квант (scaling):</strong> минимальная модель — CAS(9e,8o), классически тривиальна. Точная диагонализация (FCI) растёт комбинаторно и упирается в ~10¹⁰ детерминантов около <strong>CAS(20,20) ≈ 40 кубитов</strong> — честный кубит-порог. Фронтирный VQE — на AWS (lightning.qubit, авто-терминирование, dead-man switch как в FeMoco). Это не заявление о превосходстве, а граница применимости классики.</p>`,
      en: `<p><strong>Centre:</strong> a minimal cluster model of the radical-oxygen HAT site of Mn–Na₂WO₄/SiO₂ — an <strong>Mn=O</strong> fragment (sextet), def2-SVP, on our own validated engine (PySCF, ROHF→AVAS→CASCI). <em>Not</em> Li/MgO. Descriptor sign convention as in the "ceiling map": <code>ΔΔE‡ = barrier(C₂H₆) − barrier(CH₄)</code>, &lt;0 = over-oxidizes the C₂ product faster than methane (the BEP wall); ≳ +4 needed for 70%.</p>
           <p><strong>NOON diagnostic (basis-invariant) — the core result:</strong> at the HAT point the AVAS(Mn 3d + O 2p) natural-orbital occupations are strongly fractional: <code>[2.00, 1.70, 1.27, 1.00, 1.00, 1.00, 0.73, 0.30]</code>. The three ~1.00 values are high-spin singly-occupied orbitals (the sextet), <em>not</em> multireference; but the pairs <strong>1.70/0.30 and 1.27/0.73</strong> are genuine static correlation (≈ one–two "effectively unpaired" electrons beyond spin). Single-determinant DFT cannot reproduce this — measured, not asserted.</p>
           <p><strong>Error-bar — honest about method (failures shown, not smoothed):</strong> the qualitative conclusion is firm. A direct quantitative ΔΔE‡ correction could not be obtained, and we show why: (1) a rigid scan dragged H past the O–H equilibrium → a false TS at the scan edge (CH₄ barrier inflated, C₂H₆ went downhill — a geometry artefact); (2) AVAS picked DIFFERENT active spaces at reactant vs TS (CAS(11e,9o) vs CAS(9e,8o)) → the CASCI barrier is undefined; (3) a fixed frontier CAS(17e,12o) converges, but a small CAS loses dynamic correlation — the CASCI CH₄ barrier came out <strong>94</strong> vs <strong>11</strong> kcal/mol for DFT (an ~80 kcal/mol gap — <em>missing dynamic correlation, not selectivity physics</em>). Conclusion: a small CASCI is a <strong>diagnostic (NOON)</strong>, not a barrier; an honest quantitative correction needs a balanced multiconfigurational method (CASPT2/NEVPT2) or VQE on the real centre. That is exactly the "quantum correction" stage.</p>
           <p><strong>Where a quantum computer is needed (scaling):</strong> the minimal model is CAS(9e,8o), classically trivial. Exact diagonalization (FCI) grows combinatorially and hits ~10¹⁰ determinants around <strong>CAS(20,20) ≈ 40 qubits</strong> — the honest qubit threshold. Frontier VQE runs on AWS (lightning.qubit, auto-termination, dead-man switch as in FeMoco). Not an advantage claim — a boundary of classical applicability.</p>`
    },
    figures: [
      { src: "assets/ocm_mnw_noon.png",
        caption: { ru: "Заполнения натуральных орбиталей (NOON) в HAT-точке реального центра Mn=O: дробные значения 1.70/1.27/0.73/0.30 (красные) — настоящая статическая корреляция, которую однодетерминантный DFT описать не может. Тройка 1.00 — высокоспиновые орбитали (секстет), не мультиреференс. Реальный CASCI, минимальный кластер Mn=O, def2-SVP; базис-инвариантно.",
                   en: "Natural-orbital occupations (NOON) at the HAT point of the real Mn=O centre: the fractional values 1.70/1.27/0.73/0.30 (red) are genuine static correlation that single-determinant DFT cannot represent. The three 1.00 values are high-spin orbitals (sextet), not multireference. Real CASCI, minimal Mn=O cluster, def2-SVP; basis-invariant." } },
      { src: "assets/ocm_mnw_scaling.png",
        caption: { ru: "Где заканчивается классическая точная диагонализация и реально нужен квантовый компьютер: точный комбинаторный рост числа детерминантов с размером активного пространства; честный порог ≈ 40 кубитов. Минимальная модель сегодня классически тривиальна — это НЕ заявление о квантовом превосходстве.",
                   en: "Where classical exact diagonalization ends and a quantum computer is genuinely needed: the exact combinatorial growth of the determinant space with active-space size; honest threshold ≈ 40 qubits. The minimal model is classically trivial today — NOT a quantum-advantage claim." } }
    ]
  }
];

const GLOSSARY = [
  { term: { ru: "OCM (окислительная конденсация метана)", en: "OCM (oxidative coupling of methane)" },
    def: { ru: "Превращение метана в C₂-углеводороды (этан/этилен) с кислородом за один шаг. Упирается в ~28% потолок выхода.",
           en: "Converting methane to C₂ hydrocarbons (ethane/ethylene) with oxygen in one step. Hits a ~28% yield ceiling." } },
  { term: { ru: "Этилен (C₂H₄)", en: "Ethylene (C₂H₄)" },
    def: { ru: "Целевой продукт: «кирпичик» для полиэтилена и почти всего пластика.",
           en: "The target product: the building block of polyethylene and most plastics." } },
  { term: { ru: "Барьер (переходное состояние)", en: "Barrier (transition state)" },
    def: { ru: "Энергетическая «горка», которую надо преодолеть для реакции. Выше барьер — медленнее реакция.",
           en: "The energy 'hill' a reaction must climb. Higher barrier = slower reaction." } },
  { term: { ru: "BEP-соотношение", en: "BEP relation" },
    def: { ru: "Связь силы связи и барьера: слабее связь → ниже барьер. Поэтому продукты с более слабыми C–H переокисляются быстрее метана — источник потолка.",
           en: "Links bond strength to barrier: weaker bond → lower barrier. So products with weaker C–H bonds over-oxidize faster than methane — the source of the ceiling." } },
  { term: { ru: "Дескриптор ΔΔG‡", en: "Descriptor ΔΔG‡" },
    def: { ru: "Наш индикатор селективности = барьер(продукт) − барьер(метан). <0 — переокисляет (плохо), >0 — бережёт продукт (пробивает стену).",
           en: "Our selectivity indicator = barrier(product) − barrier(methane). <0 over-oxidizes (bad), >0 spares the product (breaks the wall)." } },
  { term: { ru: "CASCI / активное пространство", en: "CASCI / active space" },
    def: { ru: "Точный многоконфигурационный расчёт в выбранном наборе ключевых орбиталей (CAS). Нужен для разрыва связей, где обычные методы врут.",
           en: "An exact multiconfigurational calculation within a chosen set of key orbitals (CAS). Needed for bond-breaking, where ordinary methods fail." } },
  { term: { ru: "VQE", en: "VQE" },
    def: { ru: "Вариационный квантовый алгоритм поиска энергии основного состояния — мост к квантовым компьютерам для пространств, недоступных классике.",
           en: "Variational Quantum Eigensolver — a route to quantum computers for spaces beyond classical reach." } },
  { term: { ru: "Конфайнмент / молекулярное сито", en: "Confinement / molecular sieving" },
    def: { ru: "Идея отделять молекулы по размеру/форме порами. Для метана и этилена по размеру не работает — они почти одинаковы.",
           en: "Separating molecules by size/shape using pores. For methane vs. ethylene it fails by size — they're almost identical." } }
];
