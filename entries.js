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
  title:    { ru: "OCM-лаборатория", en: "OCM Lab" },
  subtitle: { ru: "Метан → этилен: открытый дневник квантово-химических экспериментов",
              en: "Methane → ethylene: an open diary of quantum-chemistry experiments" },
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
