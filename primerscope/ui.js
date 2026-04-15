// ── STATE ────────────────────────────────────────────────────────────────────
const state = {
  sequences:      [],
  consensusSeq:   '',
  seqLen:         0,
  jsdScores:      [],
  windowedScores: [],
  colDistribs:    [],
  candidates:     [],
  showFull:       false,
  activeCandidate: null,
  k:          CONFIG.k,
  c:          CONFIG.c,
  minLen:     CONFIG.minLen,
  maxLen:     CONFIG.maxLen,
  TRIM_LIMIT: CONFIG.trimLimit,
  WINDOW:     CONFIG.window,
  LAMBDA_JSD: CONFIG.lambdaJSD,
};

// ── HELPERS ──────────────────────────────────────────────────────────────────
function toast(msg, err = false) {
  const t = document.getElementById('toast');
  t.textContent = msg;
  t.className = 'toast show' + (err ? ' err' : '');
  setTimeout(() => { t.className = 'toast'; }, 2500);
}

function setStatus(msg, color = '') {
  document.getElementById('statusMsg').textContent = msg;
  document.getElementById('statusDot').className = 'status-dot' + (color ? ` ${color}` : '');
}

function syncK(v) {
  state.k = parseFloat(v);
  document.getElementById('kVal').textContent = parseFloat(v).toFixed(2);
}

function syncC(v) {
  state.c = parseFloat(v);
  document.getElementById('cVal').textContent = parseFloat(v).toFixed(2);
}

// ── ANALYSIS ─────────────────────────────────────────────────────────────────
function runAnalysis() {
  if (!state.sequences.length) { toast('No alignment loaded', true); return; }

  state.minLen = parseInt(document.getElementById('minLen').value);
  state.maxLen = parseInt(document.getElementById('maxLen').value);
  if (state.minLen < 10 || state.maxLen > 40 || state.minLen > state.maxLen) {
    toast('Invalid length range (10–40, min ≤ max)', true); return;
  }
  state.k = parseFloat(document.getElementById('kSlider').value);
  state.c = parseFloat(document.getElementById('cSlider').value);

  setStatus('Computing Window-JSD...', 'yellow');

  setTimeout(() => {
    try {
      const { raw, windowed, colDistribs } = computeJSD(state.sequences, state.WINDOW, state.LAMBDA_JSD);
      state.jsdScores      = raw;
      state.windowedScores = windowed;
      state.colDistribs    = colDistribs;

      const L = state.seqLen;
      const primers = [];
      for (let len = state.minLen; len <= state.maxLen; len++) {
        for (let start = 0; start <= L - len; start++) {
          const s = primerScore(windowed, start, len);
          if (s >= state.k) primers.push({ start, len, score: s });
        }
      }

      primers.sort((a, b) => b.score - a.score);
      const kept = [];
      for (const p of primers) {
        let overlap = false;
        for (const k2 of kept) {
          const ov = Math.min(p.start + p.len, k2.start + k2.len) - Math.max(p.start, k2.start);
          if (ov > 0.5 * Math.min(p.len, k2.len)) { overlap = true; break; }
        }
        if (!overlap) kept.push(p);
      }

      state.candidates = kept;
      renderCandidates();
      renderSeq();

      document.getElementById('statsSection').style.display = 'block';
      document.getElementById('runSummary').innerHTML =
        `Positions scored: <strong style="color:var(--accent)">${L}</strong><br>` +
        `Candidates found: <strong style="color:var(--accent)">${kept.length}</strong><br>` +
        `Score range: <strong style="color:var(--accent)">${kept.length
          ? kept[kept.length - 1].score.toFixed(3) + ' – ' + kept[0].score.toFixed(3) : '—'
        }</strong><br>` +
        `Window size: <strong style="color:var(--accent)">±${state.WINDOW} nt</strong>`;

      document.getElementById('threshDisplay').textContent = state.k.toFixed(2);
      document.getElementById('candidatesSection').style.display = 'flex';
      setStatus(`Analysis complete — ${kept.length} candidates`, 'green');
    } catch (e) {
      toast('Analysis error: ' + e.message, true);
      setStatus('Analysis error', '');
    }
  }, 10);
}

// ── DISPLAY SETTINGS ─────────────────────────────────────────────────────────
function getSeqDisplayOpts() {
  return {
    fontSize:      parseInt(document.getElementById('seqFontSlider')?.value || 14),
    colorBases:    document.getElementById('colorBasesCb')?.checked !== false,
    labelInterval: parseInt(document.getElementById('labelIntervalSel')?.value || 10),
  };
}

function applySeqDisplay() {
  const opts = getSeqDisplayOpts();
  document.getElementById('seqFontVal').textContent = opts.fontSize;
  document.getElementById('seqScroll').style.fontSize = opts.fontSize + 'px';
  renderSeq(state.activeCandidate);
}

// ── RENDERING ────────────────────────────────────────────────────────────────
function colorBase(b, useColor) {
  if (!useColor) return b;
  const cls = { A: 'base-A', T: 'base-T', G: 'base-G', C: 'base-C', '-': 'base-gap', '.': 'base-gap' }[b] || '';
  return cls ? `<span class="${cls}">${b}</span>` : b;
}

const BASES_PER_ROW = 60;

function renderSeq(highlightCand = null) {
  const seq = state.consensusSeq;
  if (!seq) return;

  const opts     = getSeqDisplayOpts();
  const limit    = state.showFull ? seq.length : Math.min(seq.length, state.TRIM_LIMIT);
  const scroll   = document.getElementById('seqScroll');
  scroll.style.fontSize   = opts.fontSize + 'px';
  scroll.style.lineHeight = '2';

  const hlSet = new Set();
  if (highlightCand) {
    for (let i = highlightCand.start; i < highlightCand.start + highlightCand.len; i++) hlSet.add(i);
  }

  const labelW = seq.length.toString().length;
  let html = '';

  for (let rowStart = 0; rowStart < limit; rowStart += BASES_PER_ROW) {
    const rowEnd = Math.min(rowStart + BASES_PER_ROW, limit);
    const rowLen = rowEnd - rowStart;

    const ruler = new Array(rowLen).fill(' ');
    for (let i = rowStart; i < rowEnd; i++) {
      const absPos = i + 1;
      if (absPos === 1 || absPos % opts.labelInterval === 0) {
        const label = String(absPos);
        const relCol = i - rowStart;
        const labelStart = relCol - label.length + 1;
        for (let j = 0; j < label.length; j++) {
          const col = labelStart + j;
          if (col >= 0 && col < rowLen) ruler[col] = label[j];
        }
      }
    }

    let basesHtml = '';
    for (let i = rowStart; i < rowEnd; i++) {
      const b     = seq[i] || '-';
      const isHL  = hlSet.has(i);
      const inner = colorBase(b, opts.colorBases);
      basesHtml += isHL
        ? `<span class="base-highlight" data-pos="${i}">${inner}</span>`
        : `<span data-pos="${i}">${inner}</span>`;
    }

    const lpad = String(rowStart + 1).padStart(labelW, ' ');
    html += `<div class="seq-row" style="align-items:baseline;">
      <span class="seq-locus" style="min-width:${labelW + 1}ch;font-size:0.8em;">${lpad}</span>
      <span class="seq-bases" style="white-space:pre;">${basesHtml}</span>
    </div>
    <div class="seq-ruler-row" style="display:flex;gap:12px;margin-bottom:2px;line-height:1;">
      <span style="min-width:${labelW + 1}ch;"></span>
      <span style="font-family:var(--mono);font-size:0.65em;color:var(--text-faint);white-space:pre;letter-spacing:0;">${ruler.join('')}</span>
    </div>`;
  }

  if (!state.showFull && seq.length > state.TRIM_LIMIT) {
    html += `<div class="seq-trimmed-note" onclick="toggleFullSeq()">▶ Showing first ${state.TRIM_LIMIT} nt of ${seq.length} — click to expand</div>`;
    const tn = document.getElementById('trimNote');
    tn.style.display = '';
    tn.textContent = '▶ Expand';
  } else if (seq.length > state.TRIM_LIMIT) {
    const tn = document.getElementById('trimNote');
    tn.textContent = '▲ Collapse';
    tn.style.display = '';
  }

  scroll.innerHTML = html;

  if (highlightCand) {
    const targetRow = Math.floor(highlightCand.start / BASES_PER_ROW);
    const rowEls    = scroll.querySelectorAll('.seq-row');
    if (rowEls[targetRow]) rowEls[targetRow].scrollIntoView({ behavior: 'smooth', block: 'center' });
  }
}

function toggleFullSeq() {
  state.showFull = !state.showFull;
  renderSeq(state.activeCandidate);
}

function renderCandidates() {
  const list    = document.getElementById('candidatesList');
  const countEl = document.getElementById('candidateCount');
  if (!state.candidates.length) {
    list.innerHTML = '<div class="no-candidates">No candidates above threshold</div>';
    countEl.textContent = '0 found';
    return;
  }
  countEl.textContent = `${state.candidates.length} found`;
  list.innerHTML = state.candidates.map((c, idx) => {
    const sc = scoreToColor(c.score);
    return `<div class="candidate-item" id="cand-${idx}" onclick="selectCandidate(${idx})">
      <span class="cand-pos">${c.start + 1}–${c.start + c.len}</span>
      <span class="cand-len">${c.len}nt</span>
      <span class="cand-seq">${getPrimerSeq(c)}</span>
      <div class="score-bar-wrap"><div class="score-bar" style="width:${(c.score * 100).toFixed(0)}%;background:${sc};"></div></div>
      <span class="cand-score" style="color:${sc};">${c.score.toFixed(3)}</span>
    </div>`;
  }).join('');
}

function scoreToColor(s) {
  if (s >= 0.85) return '#4af0a0';
  if (s >= 0.70) return '#4ab8f0';
  if (s >= 0.55) return '#f0a04a';
  return '#f04a7a';
}

function getPrimerSeq(cand) {
  return state.consensusSeq.slice(cand.start, cand.start + cand.len);
}

function selectCandidate(idx) {
  document.querySelectorAll('.candidate-item').forEach(el => el.classList.remove('active'));
  const el = document.getElementById(`cand-${idx}`);
  if (el) el.classList.add('active');
  state.activeCandidate = state.candidates[idx];
  renderSeq(state.activeCandidate);
  openDetail(idx);
}

// ── DETAIL PANEL ─────────────────────────────────────────────────────────────
function plasmaColor(t) {
  const r  = Math.round(Math.max(0, Math.min(255,
    t < 0.25 ? 13 + t*4*90 : t < 0.5 ? 103 + (t-0.25)*4*74 :
    t < 0.75 ? 177 + (t-0.5)*4*46 : 223 + (t-0.75)*4*17)));
  const g  = Math.round(Math.max(0, Math.min(255,
    t < 0.25 ? 8 + t*4*(-8) : t < 0.5 ? 0 + (t-0.25)*4*42 :
    t < 0.75 ? 42 + (t-0.5)*4*58 : 100 + (t-0.75)*4*149)));
  const b2 = Math.round(Math.max(0, Math.min(255,
    t < 0.25 ? 135 + t*4*(-5) : t < 0.5 ? 168 + (t-0.25)*4*(-24) :
    t < 0.75 ? 98 + (t-0.5)*4*(-62) : 36 + (t-0.75)*4*(-3))));
  return `rgb(${r},${g},${b2})`;
}

function openDetail(idx) {
  const cand      = state.candidates[idx];
  const primerSeq = getPrimerSeq(cand);
  const total     = state.sequences.length;

  const locusScores = [];
  for (let i = cand.start; i < cand.start + cand.len; i++) {
    locusScores.push({ pos: i + 1, jsd: state.windowedScores[i], raw: state.jsdScores[i] });
  }

  function locusVariants(pos) {
    const cnt = {};
    let ngap = 0;
    for (const s of state.sequences) {
      const b = (s.seq[pos] || '-').toUpperCase();
      if (b === '-' || b === '.') ngap++;
      else cnt[b] = (cnt[b] || 0) + 1;
    }
    return {
      cnt, nvalid: total - ngap, ngap,
      sorted: Object.entries(cnt).sort((a, b) => b[1] - a[1])
               .map(([b, n]) => ({ base: b, count: n, pct: n / total })),
    };
  }

  const aboveC   = locusScores.filter(l => l.jsd >= state.c).length;
  const cFrac    = aboveC / cand.len;
  const gcCount  = [...primerSeq].filter(b => 'GC'.includes(b)).length;
  const gcRatio  = gcCount / primerSeq.length;
  const A = [...primerSeq].filter(b => b === 'A').length;
  const T = [...primerSeq].filter(b => b === 'T').length;
  const G = [...primerSeq].filter(b => b === 'G').length;
  const C = [...primerSeq].filter(b => b === 'C').length;
  const Tm = 2*(A+T) + 4*(G+C);

  const haplotypeCounts = {};
  for (const s of state.sequences) {
    const hap = s.seq.slice(cand.start, cand.start + cand.len);
    haplotypeCounts[hap] = (haplotypeCounts[hap] || 0) + 1;
  }
  const haplotypes  = Object.entries(haplotypeCounts)
    .map(([seq, count]) => ({ seq, count, pct: count / total }))
    .sort((a, b) => b.count - a.count);
  const dominantHap = haplotypes[0]?.seq || primerSeq;

  document.getElementById('detailId').textContent =
    `PRIMER #${idx + 1}  pos ${cand.start + 1}–${cand.start + cand.len}`;

  function chip(label, value, color, warn = '') {
    return `<span style="display:inline-flex;flex-direction:column;background:var(--surface3);border:1px solid var(--border);border-radius:3px;padding:3px 9px;gap:1px;white-space:nowrap;">
      <span style="font-family:var(--mono);font-size:8px;color:var(--text-faint);letter-spacing:0.1em;">${label}</span>
      <span style="font-family:var(--mono);font-size:13px;font-weight:600;color:${color};">${value}${warn}</span>
    </span>`;
  }
  document.getElementById('detailChips').innerHTML =
    chip('SCORE',           cand.score.toFixed(4), scoreToColor(cand.score)) +
    chip('CONSERVED LOCI', `${(cFrac*100).toFixed(1)}%`, 'var(--accent)',
         `<span style="font-size:9px;color:var(--text-faint);margin-left:4px;">${aboveC}/${cand.len} ≥${state.c.toFixed(2)}</span>`) +
    chip('GC', `${(gcRatio*100).toFixed(1)}%`,
         gcRatio >= 0.4 && gcRatio <= 0.6 ? 'var(--accent)' : 'var(--accent3)',
         gcRatio >= 0.4 && gcRatio <= 0.6 ? ' ✓' : ' ⚠') +
    chip('Tm',     `${Tm}°C`,       'var(--accent2)') +
    chip('LENGTH', `${cand.len} nt`, 'var(--text)') +
    chip('DEPTH',  `${total}`,       'var(--text-dim)');

  const trackEl = document.getElementById('detailTrack');
  const rulerEl = document.getElementById('detailRuler');
  const tip     = document.getElementById('locusTooltip');
  trackEl.innerHTML = '';
  rulerEl.innerHTML = '';

  locusScores.forEach((ls, i) => {
    const base    = primerSeq[i] || '-';
    const bgColor = plasmaColor(ls.jsd);
    const textCol = ls.jsd > 0.6 ? '#111' : '#eee';
    const cell    = document.createElement('span');
    cell.className = 'locus-base-cell';
    Object.assign(cell.style, {
      display: 'inline-flex', alignItems: 'center', justifyContent: 'center',
      flex: '1', minWidth: '18px', height: '34px',
      background: bgColor, color: textCol,
      fontFamily: 'var(--mono)', fontSize: '14px', fontWeight: '700',
      borderRadius: '2px', cursor: 'crosshair',
      borderBottom: ls.jsd >= state.c ? '3px solid rgba(255,255,255,0.55)' : '3px solid transparent',
      userSelect: 'none',
    });
    cell.textContent = base;

    const vd = locusVariants(cand.start + i);
    const BASE_COLORS = { A: '#5de8a8', T: '#f07a5a', G: '#5ab4f0', C: '#e0c04a' };

    cell.addEventListener('mouseenter', () => {
      const varLines = vd.sorted.slice(0, 5).map(v =>
        `  <span style="color:${BASE_COLORS[v.base] || 'var(--text)'};">${v.base}</span>` +
        `  ${v.count.toString().padStart(5)}  <b style="color:var(--accent3)">${(v.pct*100).toFixed(2)}%</b>`
      ).join('<br>');
      const gapLine = vd.ngap > 0
        ? `<br>  <span style="color:var(--text-faint);">-</span>  ${vd.ngap.toString().padStart(5)}  <span style="color:var(--text-faint)">${(vd.ngap/total*100).toFixed(2)}%</span>`
        : '';
      tip.innerHTML =
        `<b style="color:var(--text)">Pos ${ls.pos}</b>  base: <b style="color:${BASE_COLORS[base] || 'var(--text)'};">${base}</b><br>` +
        `JSD: <b style="color:var(--accent2)">${ls.jsd.toFixed(4)}</b>  raw: ${ls.raw.toFixed(4)}<br>` +
        `<span style="color:var(--text-faint);font-size:9px;">──────────────────</span><br>` +
        varLines + gapLine;
      tip.style.display = 'block';
    });
    cell.addEventListener('mousemove', e => {
      tip.style.left = (e.clientX - tip.offsetWidth) + 'px';
      tip.style.top  = (e.clientY - tip.offsetHeight - 128) + 'px';
    });
    cell.addEventListener('mouseleave', () => { tip.style.display = 'none'; });
    trackEl.appendChild(cell);

    const relPos    = i + 1;
    const show      = (relPos === 1 || relPos % 5 === 0 || relPos === cand.len);
    const rulerSpan = document.createElement('span');
    Object.assign(rulerSpan.style, {
      display: 'inline-flex', alignItems: 'center', justifyContent: 'center',
      flex: '1', minWidth: '18px', whiteSpace: 'nowrap', overflow: 'visible',
    });
    rulerSpan.textContent = show ? ls.pos : '';
    rulerEl.appendChild(rulerSpan);
  });

  // Haplotype alignment body
  const hapCount       = haplotypes.length;
  const isMonomorphic  = hapCount === 1;
  const body           = document.getElementById('detailBody');

  let html = `<div style="padding:0 0 4px;font-family:var(--mono);font-size:9px;color:var(--text-faint);letter-spacing:0.1em;text-transform:uppercase;margin-bottom:8px;">
    Haplotype Alignment
    <span style="color:var(--text-faint);font-weight:400;margin-left:8px;font-size:9px;">
      ${hapCount} unique haplotype${hapCount !== 1 ? 's' : ''} · ${total} sequences · BLAST-style
    </span>
  </div>`;

  if (isMonomorphic) {
    html += `<div style="color:var(--text-faint);font-family:var(--mono);font-size:11px;margin-bottom:16px;">
      All ${total} sequences identical at this primer window
    </div>`;
  } else {
    const domSeq   = haplotypes[0].seq;
    const top3     = haplotypes.slice(0, 3);
    const restHaps = haplotypes.slice(3);
    const moreId   = 'hap-more-' + cand.start;

    function blastAlign(qSeq, sSeq, sRank, sPct, sCount) {
      let identLine = '', colorQ = '', colorS = '';
      const BASE_COLORS_LOCAL = { A: '#5de8a8', T: '#f07a5a', G: '#5ab4f0', C: '#e0c04a' };
      for (let i = 0; i < Math.max(qSeq.length, sSeq.length); i++) {
        const q = qSeq[i] || '-', s = sSeq[i] || '-';
        const same = q === s, gap = q === '-' || s === '-';
        identLine += same ? '|' : (gap ? ' ' : '.');
        function bspan(b, mis) {
          if (b === '-') return '<span style="color:var(--text-faint)">-</span>';
          const bc = BASE_COLORS_LOCAL[b] || 'var(--text)';
          return mis
            ? `<span style="color:var(--accent3);background:rgba(240,160,74,0.2);font-weight:700;border-radius:1px">${b}</span>`
            : `<span style="color:${bc}">${b}</span>`;
        }
        colorQ += bspan(q, !same && !gap);
        colorS += bspan(s, !same && !gap);
      }
      const nId  = [...qSeq].filter((c, i) => c === sSeq[i]).length;
      const pId  = (nId / qSeq.length * 100).toFixed(1);
      const nMis = qSeq.length - nId;
      return `<div style="background:var(--surface3);border:1px solid var(--border);border-radius:3px;padding:10px 14px;margin-bottom:8px;font-family:var(--mono);font-size:12px;">
        <div style="display:flex;justify-content:space-between;margin-bottom:7px;">
          <span style="color:var(--accent3);font-size:10px;font-weight:600;">Haplotype #${sRank}</span>
          <span style="color:var(--text-dim);font-size:10px;">${sCount} seq · ${sPct}% · <span style="color:var(--accent)">Identity: ${pId}%</span> · <span style="color:var(--accent4)">Mismatches: ${nMis}</span></span>
        </div>
        <div style="overflow-x:auto;"><div style="display:inline-block;white-space:pre;line-height:1.7;">` +
        `<span style="color:var(--text-faint);font-size:10px;display:inline-block;width:64px;text-align:right;padding-right:8px;">Query  1</span><span style="letter-spacing:0.06em;">${colorQ}</span><span style="color:var(--text-faint);font-size:10px;padding-left:6px;">${qSeq.length}</span>\n` +
        `<span style="display:inline-block;width:64px;"></span><span style="color:var(--text-faint);letter-spacing:0.06em;">${identLine}</span>\n` +
        `<span style="color:var(--text-faint);font-size:10px;display:inline-block;width:64px;text-align:right;padding-right:8px;">Sbjct  1</span><span style="letter-spacing:0.06em;">${colorS}</span><span style="color:var(--text-faint);font-size:10px;padding-left:6px;">${sSeq.length}</span>` +
        `</div></div></div>`;
    }

    const domColorSeq = [...domSeq].map(b => {
      const c = { A: '#5de8a8', T: '#f07a5a', G: '#5ab4f0', C: '#e0c04a' }[b] || 'var(--text)';
      return `<span style="color:${c}">${b}</span>`;
    }).join('');

    html += `<div style="background:rgba(74,240,160,0.05);border:1px solid rgba(74,240,160,0.2);border-radius:3px;padding:10px 14px;margin-bottom:10px;font-family:var(--mono);font-size:12px;">
      <div style="display:flex;justify-content:space-between;margin-bottom:6px;">
        <span style="color:var(--accent);font-size:10px;font-weight:600;">Dominant Haplotype #1 (Query)</span>
        <span style="color:var(--text-dim);font-size:10px;">${haplotypes[0].count} seq · ${(haplotypes[0].pct*100).toFixed(1)}%</span>
      </div>
      <div style="overflow-x:auto;white-space:pre;">
        <span style="color:var(--text-faint);font-size:10px;display:inline-block;width:64px;text-align:right;padding-right:8px;">Query  1</span>
        <span style="letter-spacing:0.06em;">${domColorSeq}</span>
        <span style="color:var(--text-faint);font-size:10px;padding-left:6px;">${domSeq.length}</span>
      </div>
    </div>`;

    top3.slice(1).forEach((h, i) => {
      html += blastAlign(domSeq, h.seq, i + 2, (h.pct*100).toFixed(1), h.count);
    });

    if (restHaps.length > 0) {
      html += `<button class="expand-variants-btn" onclick="toggleVars('${moreId}')">▸ Show ${restHaps.length} more haplotype${restHaps.length > 1 ? 's' : ''}</button>
      <div id="${moreId}" style="display:none;margin-top:8px;">`;
      restHaps.forEach((h, i) => {
        html += blastAlign(domSeq, h.seq, i + 4, (h.pct*100).toFixed(1), h.count);
      });
      html += `</div>`;
    }
  }

  body.innerHTML = html;
  document.getElementById('detailPanel').classList.add('open');
}

function toggleVars(id) {
  const el = document.getElementById(id);
  if (el) el.style.display = el.style.display === 'none' ? 'block' : 'none';
}

function closeDetail() {
  document.getElementById('detailPanel').classList.remove('open');
  document.querySelectorAll('.candidate-item').forEach(el => el.classList.remove('active'));
  state.activeCandidate = null;
  renderSeq();
}

// ── FILE LOAD ────────────────────────────────────────────────────────────────
function handleFile(file) {
  const reader = new FileReader();
  reader.onload = e => {
    const seqs = parseFASTA(e.target.result);
    if (!seqs.length) { toast('No valid FASTA sequences found', true); return; }

    const lens = [...new Set(seqs.map(s => s.seq.length))];
    let seqsAligned = seqs;
    if (lens.length > 1) {
      toast(`Warning: unequal lengths. Using first ${lens[0]} nt.`);
      seqsAligned = seqs.map(s => ({ ...s, seq: s.seq.slice(0, lens[0]) }));
    }

    state.sequences      = seqsAligned;
    state.seqLen         = seqsAligned[0].seq.length;
    state.consensusSeq   = buildConsensus(seqsAligned);
    state.candidates     = [];
    state.activeCandidate = null;
    state.showFull       = false;

    document.getElementById('uploadZone').style.display    = 'none';
    document.getElementById('seqDisplayWrap').style.display = 'flex';
    document.getElementById('nSeqDisplay').textContent     = seqsAligned.length;
    document.getElementById('lenDisplay').textContent      = `${state.seqLen} nt`;
    document.getElementById('seqLabel').className          = 'panel-label active';
    document.getElementById('runBtn').disabled             = false;
    document.getElementById('candidatesSection').style.display = 'none';
    document.getElementById('statsSection').style.display      = 'none';

    renderSeq();
    setStatus(`Loaded ${seqsAligned.length} sequences × ${state.seqLen} nt`, 'blue');
    toast(`Loaded: ${seqsAligned.length} sequences`);
  };
  reader.readAsText(file);
}

// ── EVENT BINDINGS ────────────────────────────────────────────────────────────
const _fileInput  = document.getElementById('fileInput');
const _uploadZone = document.getElementById('uploadZone');

_fileInput.addEventListener('change', e => {
  if (e.target.files[0]) handleFile(e.target.files[0]);
});

_uploadZone.addEventListener('dragover', e => {
  e.preventDefault(); e.stopPropagation();
  _uploadZone.classList.add('drag-over');
});
_uploadZone.addEventListener('dragleave', e => {
  e.stopPropagation();
  _uploadZone.classList.remove('drag-over');
});
_uploadZone.addEventListener('drop', e => {
  e.preventDefault(); e.stopPropagation();
  _uploadZone.classList.remove('drag-over');
  if (e.dataTransfer.files[0]) handleFile(e.dataTransfer.files[0]);
});

// Init sliders from CONFIG
syncK(document.getElementById('kSlider').value);
syncC(document.getElementById('cSlider').value);
