# 819-Day Station: Palenque Temple of the Sun → Temple of the Inscriptions

Session log — verification notes and numerical explorations tied to the 819-day count text on the left side of the Temple of the Sun's central tablet at Palenque.

## The two dates

**Base 819-day station:** Long Count 1.6.14.11.2, 1 Ik 10 Tzec
- GMT correlation (584283): JD 776745
- Julian calendar: August 12, 2587 BCE
- Proleptic Gregorian: July 22, 2587 BCE

**Medallion Series date, Temple of the Inscriptions, 63 Calendar Rounds later:** Long Count 9.12.16.2.2
- GMT correlation (584283): JD 1972485
- Julian calendar: May 15, 688 CE

## Long Count math (verified)

Using GMT correlation constant 584283:
- 1.6.14.11.2 → JD 776745
- 9.12.16.2.2 → JD 1972485
- Difference: 1,195,740 days = exactly 63 × 18,980 (Calendar Round length) — confirms "sixty-three Calendar Rounds later" exactly.
- JD 1972485 in the *Julian* calendar (the convention Mayanist literature uses for pre-1582 dates) falls on May 15, 688 — matches the cited date exactly.

## Astronomical verification for 5/15/688 (Julian calendar), computed via VSOP87 geocentric ecliptic longitudes

- **Jupiter–Mars separation:** 5.2° — not exact conjunction, but close and widening. The actual conjunction (0.9° separation) occurred ~31 days earlier, April 14, 688.
- **Jupiter & Mars motion:** both direct (forward) on this date. Jupiter's prior retrograde-to-direct station was ~Jan 2, 688; Mars's was ~Jan 30, 688 — so both had resumed forward motion well before 5/15/688, consistent with (if not simultaneous with) "beginning forward motion after retrograde."
- **Saturn:** retrograde on 5/15/688 (~−0.22°/day), approaching its second stationary point (retrograde→direct) computed at JD 1972523.8 → **June 23, 688 (Julian calendar)** — an exact match to the source text's "(6/23/688; Meeus n.d.)" citation.

## Context

Per the source text: Lounsbury did not note the Saturn retrograde/station detail. The hypothesis flagged is that Jupiter and Saturn's synodic cycles are linked through the 819-day count.

## Jupiter's revolution: sidereal vs. synodic (two different numbers)

Confirmed there are genuinely two numbers, not one:
- **Sidereal period ("star to star"):** true orbital period against the fixed stars. ≈ **11.862 years ≈ 4332.59 days**. Not directly observable in one sitting — requires tracking Jupiter's position against the zodiac over years and correcting for Earth's own motion.
- **Synodic period ("horizon to horizon" — conjunction to conjunction with the Sun, heliacal rising to heliacal rising, station to station, opposition to opposition):** ≈ **398.88 days ≈ 1.092 years**. This is the number a naked-eye observer can actually count directly, without needing a multi-year correction.
- Related by 1/T_syn = 1/T_earth − 1/T_jupiter. Quick VSOP87 check of individual real Jupiter–Sun conjunctions in one stretch gave ~399–403 days for single cycles (vs. the 398.88-day mean) — normal cycle-to-cycle variation from orbital eccentricity.
- For a naked-eye/day-count tradition, the synodic period (~399 d) is the directly countable quantity; the sidereal period (4332.59 d) is a derived, longer-baseline quantity.

## 819-day count: numerical commensurations with planetary periods

Systematic small-integer search (n × 819 days ≈ m × planetary period), cross-validated against the already-published Mars relation before trusting it on Jupiter/Saturn:

- **Mars (known/published relation):** 20 × 819 = 16,380 d ≈ 21 × 779.94-d Mars synodic period (16,378.7 d) — error **+0.008%**. Matches literature (same 780-day round used in the Dresden Codex Mars table).
- **Jupiter:** 19 × 819 = 15,561 d ≈ 39 × 398.88-d Jupiter synodic period (15,556.3 d) — error **+0.030%**.
- **Saturn:** 6 × 819 = 4,914 d ≈ 13 × 378.09-d Saturn synodic period (4,915.2 d) — error **−0.024%**. Scales cleanly (12×819≈26×Sat-syn, 18×819≈39×Sat-syn, same tightness).
- **Jupiter + Saturn jointly, same multiple:** 18 × 819 = 14,742 d (40.36 yr) ≈ 37 Jupiter synodic periods (−0.11%) *and* ≈ 39 Saturn synodic periods (−0.02%) simultaneously.
- **Tightest relation found — Jupiter–Saturn great-conjunction cycle:** 62 × 819 = 50,778 d ≈ **7 × 7,253.46-d Jupiter–Saturn synodic (great-conjunction) period** (50,774.2 d) — error **+0.007%**. That's 139.02 years = 7 Trigon cycles, landing almost exactly on 62 whole 819-day stations. Tighter than the Mars relation.

Caveat: this is a computational small-integer search (numerology in the neutral sense), not evidence of intent — but the Mars case cross-validates the method against an already-accepted result, and the Jupiter–Saturn great-conjunction match is the tightest of all found.

## Jupiter's revolution: user's count-derived formula vs. tropical/sidereal correction

User-supplied formula: 63 × 260 × 146 / 552, proposed as a "solar-centered" (short) estimate of Jupiter's revolution, vs. the standard "star to star" sidereal figure (4332.59 d), with the gap attributed to a tropical/sidereal correction factor of 25920/25919 (25,920-year precession cycle).

Verified arithmetic (exact):
- 63 × 260 × 146 / 552 = 99645/23 = **4332.391304... days**
- × 25920/25919 = **4332.558456 days**
- Difference from commonly cited 4332.59: **−0.0315 days (0.0007%)**

Notes/caveats logged:
- The tropical-vs-sidereal distinction is real: a frame that tracks the moving equinox vs. one fixed against the stars will disagree, and the gap compounds over time. A quick independent VSOP87 check (Jupiter's heliocentric longitude "of date" across one orbit) did come back shorter than 4332.59 days, qualitatively consistent — but single-orbit sampling wasn't precise enough to pin the exact figure, so treat that check as suggestive, not confirmatory.
- On the 25,920-year figure specifically: the **modern measured axial precession period is ≈25,772 years**, not 25,920. Back-solving the precession period implied by the standard tropical year (365.24219 d) and sidereal year (365.256363 d) gives ≈25,770 years — matching the modern measurement, not 25,920. The 25,920 figure is the traditional "Great Year"/Platonic Year approximation (360° ÷ 72 yr/degree), attractive because it divides evenly by 12/30/60/360, but it's ~150 years short of the instrumentally measured (VLBI, lunar laser ranging) modern value. This is a normal historical refinement (same as the tropical year, obliquity, and lunar month getting refined digits over time), not evidence of suppression.
- At this precision, the observed 0.03-day gap between the formula and 4332.59 is smaller than the difference the 25920 vs. ~25770 precession constants would themselves produce — so the match doesn't actually discriminate between the two values; it's consistent with either.

## Open thread

Ongoing exploration of what else the 819-day station encodes — to be extended as the discussion continues.
