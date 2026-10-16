// Builds "Sirius - Proper Motion, Companion and Calendar Cycle.docx" from the charts made in this folder.
const fs = require('fs'), path = require('path');
const { Document, Packer, Paragraph, TextRun, ImageRun, HeadingLevel, AlignmentType, Table, TableRow, TableCell,
        WidthType, ShadingType, BorderStyle, LevelFormat, PageBreak, Footer, PageNumber } = require('docx');
const here = __dirname;
const OUT = path.join(here, 'Sirius - Proper Motion, Companion and Calendar Cycle.docx');
const W = 9360; // text width in DXA (US Letter, 1-inch margins)

const P = (text, opts = {}) => new Paragraph({ spacing: { after: 140, line: 288 }, ...opts,
  children: (Array.isArray(text) ? text : [text]).map(t => typeof t === 'string' ? new TextRun(t) : t) });
const B = t => new TextRun({ text: t, bold: true });
const I = t => new TextRun({ text: t, italics: true });
const H1 = t => new Paragraph({ heading: HeadingLevel.HEADING_1, spacing: { before: 320, after: 160 }, children: [new TextRun(t)] });
const H2 = t => new Paragraph({ heading: HeadingLevel.HEADING_2, spacing: { before: 240, after: 120 }, children: [new TextRun(t)] });
const bullet = (parts) => new Paragraph({ numbering: { reference: 'bul', level: 0 }, spacing: { after: 80, line: 276 },
  children: (Array.isArray(parts) ? parts : [parts]).map(t => typeof t === 'string' ? new TextRun(t) : t) });
function fig(file, wpx, hpx, caption, w = 624) {
  const h = Math.round(w * hpx / wpx); // 6.5 in at 96 dpi
  return [new Paragraph({ alignment: AlignmentType.CENTER, spacing: { before: 120, after: 60 }, keepNext: true,
            children: [new ImageRun({ type: 'png', data: fs.readFileSync(path.join(here, file)), transformation: { width: w, height: h },
                                      altText: { title: caption, description: caption, name: file } })] }),
          new Paragraph({ alignment: AlignmentType.CENTER, spacing: { after: 220 },
            children: [new TextRun({ text: caption, italics: true, size: 18, color: '555555' })] })];
}
const border = { style: BorderStyle.SINGLE, size: 4, color: 'BFC7D5' };
const borders = { top: border, bottom: border, left: border, right: border };
function table(rows, widths) {
  return new Table({ width: { size: W, type: WidthType.DXA }, columnWidths: widths,
    rows: rows.map((r, ri) => new TableRow({ tableHeader: ri === 0, children: r.map((c, ci) => new TableCell({
      borders, width: { size: widths[ci], type: WidthType.DXA }, margins: { top: 60, bottom: 60, left: 100, right: 100 },
      shading: ri === 0 ? { fill: '1F3A5F', type: ShadingType.CLEAR, color: 'auto' } : (ri % 2 === 0 ? { fill: 'F2F5F9', type: ShadingType.CLEAR, color: 'auto' } : undefined),
      children: [new Paragraph({ alignment: ci === 0 ? AlignmentType.LEFT : AlignmentType.CENTER,
        children: [new TextRun({ text: c, bold: ri === 0, color: ri === 0 ? 'FFFFFF' : '222222', size: 19 })] })] })) })) });
}
const spacer = () => new Paragraph({ spacing: { after: 120 }, children: [] });

const body = [
  new Paragraph({ alignment: AlignmentType.CENTER, spacing: { after: 80 }, children: [new TextRun({ text: 'Sirius', bold: true, size: 52, color: '1F3A5F' })] }),
  new Paragraph({ alignment: AlignmentType.CENTER, spacing: { after: 80 }, children: [new TextRun({ text: 'Proper Motion, Companion and Calendar Cycle', size: 30, color: '1F3A5F' })] }),
  new Paragraph({ alignment: AlignmentType.CENTER, spacing: { after: 360 }, children: [new TextRun({ text: 'A companion note to the Tzolkin clock · October 2026', size: 20, color: '777777' })] }),

  H1('1. Why Sirius'),
  P('Sirius is the brightest star in the night sky and one of the nearest, 8.6 light years away. It has been a timekeeper for at least five thousand years: in Egypt its first appearance in the dawn twilight, the heliacal rising, announced the Nile flood and was tied to the start of the year. This note gathers what the charts in the Stars folder show about three questions:'),
  bullet('How does Sirius itself move across the sky, over 500 and over 100,000 years?'),
  bullet('How does its white-dwarf companion, Sirius B, pull it off a straight path?'),
  bullet('Is the yearly cycle of Sirius still 365.25 days, and when was that figure exact?'),
  P('Procyon, the other bright nearby star with a white-dwarf companion, is compared along the way.'),

  H2('Sirius at a glance'),
  table([
    ['Quantity', 'Value'],
    ['Distance', '8.60 light years (parallax 379.21 mas, Hipparcos 2007)'],
    ['Brightness', 'magnitude −1.46 (Sirius A, type A1 V)'],
    ['Masses', 'Sirius A 2.06 suns; Sirius B 1.02 suns (white dwarf)'],
    ['Proper motion', '1.339″ a year toward position angle 204° (south-southwest)'],
    ['Radial velocity', '−5.5 km/s (approaching)'],
    ['Speed relative to the Sun', '17.6 km/s'],
    ['Orbit of Sirius B', '50.13 years; 7.50″ (19.8 AU) semi-major axis; eccentricity 0.59'],
  ], [3000, 6360]),
  spacer(),

  H1('2. Proper motion: Sirius on its own'),
  P(['A star’s ', I('proper motion'), ' is its own slow drift across the sky, caused by its real motion through space relative to the Sun. It is separate from precession of the equinoxes, which shifts the whole coordinate grid and moves every star’s coordinates together. All of the sky charts here use fixed J2000 coordinates, so precession is left out and only the star’s own motion shows.']),
  P(['Sirius moves 1.339 arcseconds a year, one of the larger proper motions among bright stars. Edmond Halley noticed in 1718 that Sirius, Arcturus and Aldebaran had shifted since Ptolemy’s catalogue; these were the first proper motions ever detected.']),
  H2('Over 100,000 years'),
  P('Run back in a straight line, Sirius was 28° from its present place 100,000 years ago, near right ascension 112.5° and declination +9.2°. It was then 11.8 light years away and a little fainter, about magnitude −0.77. Procyon also moves quickly, about 35° in the same span.'),
  ...fig('sirius_motion_100000_years.png', 2248, 1153, 'Figure 1. Sirius’s path over the last 100,000 years, with its distance and brightness seen from Earth.'),
  H2('Over 500 years'),
  P('Zoomed in to the last five centuries, Sirius moved 11.2 arcminutes from 1526 to 2026, about a third of the width of the full Moon. One arcminute takes about 45 years; each ten-year step is about 13 arcseconds.'),
  ...fig('sirius_motion_500_years_zoom.png', 1342, 1446, 'Figure 2. Sirius’s proper motion from 1526 to 2026, with ten-year steps. East is to the left, as on a sky map.', 430),

  new Paragraph({ children: [new PageBreak()] }),
  H1('3. The pull of Sirius B'),
  P(['Sirius is a double star. What moves in a perfectly straight line is not Sirius A but the ', B('centre of mass'), ' of the pair. Sirius A circles that point once every 50.13 years, opposite Sirius B. Because Sirius B has half Sirius A’s mass, Sirius A swings about a third as far as Sirius B does.']),
  bullet([B('Size of the wobble: '), 'Sirius A is never more than 3.7″ from the centre of mass. Along its direction of travel it runs up to 3.1″ ahead and 1.2″ behind; to the side it strays about 3.6″ from end to end.']),
  bullet([B('Discovery: '), 'in 1844 Friedrich Bessel found, from a century of positions, that Sirius speeds up and slows down along its track, and concluded that an unseen companion was pulling it. The two stars were then near their closest approach.']),
  bullet([B('First sight: '), 'in 1862 Alvan Graham Clark saw Sirius B while testing a new 18½-inch telescope lens, about 10″ from Sirius A.']),
  bullet([B('Closest and farthest: '), 'the stars come within 8.1 AU of each other (1844, 1894, 1994) and separate to 31.5 AU (1869, 1919, 2019).']),
  P('On the 500-year chart the wobble is only a ripple, but the close-up of one orbit and the panel of offsets by year show it clearly.'),
  ...fig('sirius_motion_500_years_with_sirius_b.png', 2510, 1533, 'Figure 3. Sirius A’s path including the pull of Sirius B: the whole 500 years (left), one orbit close up, both orbits about the centre of mass, and the wobble year by year.'),

  H2('How long is Sirius B’s orbit? True and “synodic” periods'),
  P(['The ', B('true orbital period'), ' is 50.128 ± 0.004 years (about 18,309 days; Bond et al. 2017), timed from one periastron (closest approach) to the next, for example 1994.57 to 2044.70. Timed this way the period does not depend on the coordinate frame or on where it is observed from.']),
  P('A binary star has no synodic period in the way a planet does, because the Earth is far outside the system. Three small effects make the period seen from here differ slightly from the true one:'),
  table([
    ['Measured against', 'Period', 'Why it differs'],
    ['Sirius’s own frame (the true clock)', '50.1284 years', '—'],
    ['As seen from Earth (light-time)', 'about 8 hours shorter', 'Sirius approaches at 5.5 km/s, so each orbit’s light arrives slightly sooner; the period is compressed by 0.0018%.'],
    ['Relative to the Sun–Sirius line', 'within about ±1 day', 'Sirius’s proper motion turns our line of sight by 1.34″ a year, a 0.005% effect.'],
    ['Back to the same position angle in coordinates of date', 'about 50.17 years (≈ 15 days longer)', 'Precession turns the coordinate grid by about 20.5″ a year at Sirius. Sirius B moves clockwise (position angle decreasing), so it takes extra time to catch up to the same angle.'],
  ], [2900, 2000, 4460]),
  spacer(),
  P(['The 50.17-year figure is the one that behaves like a synodic period: the return of Sirius B to the same position angle in a precessing (equinox-of-date) frame, just as a planet’s synodic period is timed against the moving Sun. In fixed J2000 coordinates, as used on the charts, that 15-day difference disappears. The published 50.13-year period is therefore already what is observed from Earth; the 8-hour light-time difference is far smaller than its ±1.6-day uncertainty.']),

  H1('4. Seen from the north galactic pole'),
  P('Looking down on the Milky Way’s disc from the north, the same motion appears as a path through space. Sirius moves 17.6 km/s relative to the Sun: 13.1 km/s along the galactic plane and 11.8 km/s southward, below the plane. In each 50-year orbit the pair advances about 140 AU across the plane, while the two stars circle 8 to 31.5 AU apart. From the galactic pole the orbit is seen almost edge-on, as a tall, thin ellipse.'),
  ...fig('sirius_motion_from_north_galactic_pole.png', 2441, 1359, 'Figure 4. Sirius A and B seen from the north galactic pole: the path from 1830 to 1880, both orbits with the steady motion removed, and the whole 500 years. Distances in AU.'),
  P([I('Note: '), I('from Earth’s sky alone one cannot tell which side of the tilted orbit faces us. The usual radial-velocity convention is used; if it were reversed, the middle panel would be mirrored, but its sizes would not change.')]),

  new Paragraph({ children: [new PageBreak()] }),
  H1('5. Sirius and Procyon compared'),
  P('Procyon, 11.5 light years away, is the only other bright neighbour with a white-dwarf companion. Bessel inferred both companions in the same 1844 paper. Arthur Auwers computed Procyon’s orbit in 1862, and John Schaeberle first saw Procyon B in 1896 at Lick Observatory.'),
  table([
    ['', 'Sirius', 'Procyon'],
    ['Distance', '8.60 ly', '11.46 ly'],
    ['Masses (suns)', '2.06 + 1.02', '1.48 + 0.59'],
    ['Orbit', '50.13 years', '40.84 years'],
    ['Eccentricity', '0.59', '0.40'],
    ['Closest – farthest', '8.1 – 31.5 AU', '9.1 – 21.2 AU'],
    ['Tilt of orbit from face-on', '44°', '31°'],
    ['Proper motion', '1.339″ a year', '1.259″ a year'],
    ['Moved in 500 years', '11.2′', '10.5′'],
    ['Wobble of the bright star', 'up to 3.7″', 'up to 1.5″'],
    ['Speed relative to the Sun', '17.6 km/s', '21.2 km/s'],
    ['Travel across the plane per orbit', '139 AU', '87 AU'],
  ], [3560, 2900, 2900]),
  spacer(),
  P('Their proper motions are almost the same, but Sirius A wobbles about two and a half times as much, because Sirius B is heavier and its orbit wider and more eccentric. Procyon is the faster star, yet more of its speed carries it south out of the galactic plane, so it covers less ground across the plane in each orbit.'),
  ...fig('sirius_vs_procyon.png', 2506, 1591, 'Figure 5. Sirius and Procyon at matching scales: their orbits on the sky (arcseconds) and from the galactic pole (AU), and how far each bright star runs ahead of or behind its steady motion.'),

  new Paragraph({ children: [new PageBreak()] }),
  H1('6. The yearly cycle of Sirius: is it still 365.25 days?'),
  P('“The year of Sirius” can mean three different things, and each has its own length:'),
  bullet([B('Heliacal rising: '), 'the first morning Sirius is seen in the dawn after its weeks hidden in the Sun’s glare. This is the Egyptian “Sothic” year, and it depends on the observer’s latitude.']),
  bullet([B('Midnight transit: '), 'Sirius crossing the meridian at the same clock time each year. Today it culminates at midnight around 1 January.']),
  bullet([B('Conjunction: '), 'the Sun passing Sirius in ecliptic longitude. This is essentially the sidereal year.']),
  P('The table gives each cycle’s length, averaged over ±100 years, from a model that includes precession, Sirius’s own motion and the slowing of the Earth’s rotation.'),
  table([
    ['Cycle (days)', '4000 BC', '2600 BC', 'AD 139', 'Today'],
    ['Heliacal rising at Memphis', '365.2503', '365.2505', '365.2511', '365.2519'],
    ['Midnight meridian transit', '365.2544', '365.2537', '365.2529', '365.2527'],
    ['Sun passes Sirius', '365.2566', '365.2564', '365.2562', '365.2560'],
    ['Vernal-equinox year (reference)', '365.2423', '365.2422', '365.2423', '365.2424'],
  ], [3360, 1500, 1500, 1500, 1500]),
  spacer(),
  H2('When 365.25 days was right'),
  P(['The heliacal rising year at Memphis stayed within about a minute of 365.25 days from before 5000 BC until about 1300 BC. That is why the Egyptian ', B('Sothic cycle'), ' works: a civil calendar of exactly 365 days slips one day every four years against Sirius, and comes back into step after 1,460 Julian years (1,461 Egyptian years). The Roman writer Censorinus recorded that the rising fell on the Egyptian New Year, 1 Thoth, in AD 139; the model puts it on 19 July (Julian) that year at Memphis, matching his record.']),
  P('In the Julian calendar, which counts exactly 365.25 days a year, the rising date has moved only from about 16 July to 22 July between 5000 BC and today. Since antiquity the rising year has slowly lengthened and is now about 2.7 minutes longer than 365.25 days.'),
  H2('What it means for the clock'),
  bullet('The midnight transit of Sirius runs 365.2527 days, about 4 minutes longer than 365.25 days.'),
  bullet('Against the Gregorian calendar (365.2425 days) the midnight transit drifts about one day later every 100 years.'),
  bullet('Against a 365.25-day count it drifts only about one day every 370 years, so 365.25 remains a good whole-number rule for Sirius.'),
  ...fig('sirius_synodic_cycle_over_time.png', 2211, 1541, 'Figure 6. Length of Sirius’s yearly cycles from 5000 BC to AD 3000 (top), and the Julian-calendar date of its first dawn appearance at Memphis (bottom).'),

  H1('7. How the figures were made'),
  bullet([B('Positions and motions: '), 'Hipparcos new reduction (van Leeuwen 2007) for position, parallax and proper motion; published radial velocities (Sirius −5.5 km/s, Procyon −3.2 km/s). Motion is a straight line in three dimensions.']),
  bullet([B('Orbits: '), 'Sirius B from Bond et al. (2017); Procyon B from Bond et al. (2015). The Hipparcos proper motions are taken as the steady motion of each pair’s centre of mass, since Hipparcos fitted these binaries with their known orbits.']),
  bullet([B('Precession: '), 'the long-term model of Vondrák et al. (2011); the Sun from Meeus’s solar theory; the Earth’s rotation from the Morrison & Stephenson ΔT formula.']),
  bullet([B('Heliacal rising: '), 'Sirius 0.5° above the horizon while the Sun is 9° below, at Memphis (29.85° N). Changing this limit, or moving to Alexandria or Thebes, changes the cycle lengths by less than 0.0004 days.']),
  bullet([B('Scripts: '), 'make_sirius_motion.py, make_sirius_500.py, make_sirius_500_zoom.py, make_sirius_500_with_b.py, make_sirius_ngp_view.py, binary_stars.py, make_binary_charts.py, make_sirius_vs_procyon.py, sothic_model.py, make_sothic_series.py and make_sothic_chart.py, all in the Stars folder.']),
];

const doc = new Document({
  creator: 'Joshua Harmon', title: 'Sirius - Proper Motion, Companion and Calendar Cycle',
  styles: {
    default: { document: { run: { font: 'Calibri', size: 22 } } },
    paragraphStyles: [
      { id: 'Heading1', name: 'Heading 1', basedOn: 'Normal', next: 'Normal', quickFormat: true, run: { size: 32, bold: true, color: '1F3A5F' }, paragraph: { outlineLevel: 0 } },
      { id: 'Heading2', name: 'Heading 2', basedOn: 'Normal', next: 'Normal', quickFormat: true, run: { size: 26, bold: true, color: '2E5A88' }, paragraph: { outlineLevel: 1 } },
    ],
  },
  numbering: { config: [{ reference: 'bul', levels: [{ level: 0, format: LevelFormat.BULLET, text: '•', alignment: AlignmentType.LEFT,
    style: { paragraph: { indent: { left: 540, hanging: 270 } } } }] }] },
  sections: [{
    properties: { page: { size: { width: 12240, height: 15840 }, margin: { top: 1440, bottom: 1440, left: 1440, right: 1440 } } },
    footers: { default: new Footer({ children: [new Paragraph({ alignment: AlignmentType.CENTER,
      children: [new TextRun({ text: 'Sirius · page ', size: 16, color: '888888' }), new TextRun({ children: [PageNumber.CURRENT], size: 16, color: '888888' })] })] }) },
    children: body,
  }],
});
Packer.toBuffer(doc).then(b => { fs.writeFileSync(OUT, b); console.log('wrote', OUT, b.length); });
