// ── STATE ────────────────────────────────────────────────────────────────────
const state = {
  sequences:      [],
  consensusSeq:   '',
  seqLen:         0,
  jsdScores:      [],
  windowedScores: [],
  colDistribs:    [],
  dominantStats:  [],
  candidates:     [],
  showFull:       false,
  activeCandidate: null,
  orientationStats: { flipped: [] },
  dominantZoom: { start: 1, end: null, barWidth: null },
  hoveredCandidate: null,
  k:          CONFIG.k,
  c:          CONFIG.c,
  dominantThreshold: CONFIG.dominantThreshold,
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

const DOMINANT_THRESHOLD_MIN = 0.1;
const DOMINANT_THRESHOLD_MAX = 0.9999999;
const DOMINANT_THRESHOLD_MAX_GAP = 1 - DOMINANT_THRESHOLD_MAX;
const DOMINANT_THRESHOLD_MIN_GAP = 1 - DOMINANT_THRESHOLD_MIN;

function sliderToDominantThreshold(v) {
  const x = Math.max(0, Math.min(1, parseFloat(v)));
  const gap = DOMINANT_THRESHOLD_MIN_GAP *
    Math.pow(DOMINANT_THRESHOLD_MAX_GAP / DOMINANT_THRESHOLD_MIN_GAP, x);
  return 1 - gap;
}

function dominantThresholdToSlider(threshold) {
  const t = Math.max(DOMINANT_THRESHOLD_MIN, Math.min(DOMINANT_THRESHOLD_MAX, parseFloat(threshold)));
  const gap = 1 - t;
  return Math.max(0, Math.min(1,
    Math.log(gap / DOMINANT_THRESHOLD_MIN_GAP) /
    Math.log(DOMINANT_THRESHOLD_MAX_GAP / DOMINANT_THRESHOLD_MIN_GAP)
  ));
}

function formatDominantThreshold(value) {
  const t = Math.max(DOMINANT_THRESHOLD_MIN, Math.min(DOMINANT_THRESHOLD_MAX, value));
  const gap = 1 - t;
  if (gap < 0.001) return `1-${gap.toExponential(2)}`;
  return t.toFixed(3);
}

function parseDominantThresholdInput(raw) {
  const text = String(raw).trim().toLowerCase();
  if (!text) return NaN;
  const oneMinus = text.match(/^1\s*-\s*(.+)$/);
  if (oneMinus) return 1 - Number(oneMinus[1]);
  return Number(text);
}

function syncDominantThreshold(v) {
  state.dominantThreshold = sliderToDominantThreshold(v);
  document.getElementById('dominantThresholdVal').value = formatDominantThreshold(state.dominantThreshold);
  renderDominantHistogram();
}

function syncDominantThresholdInput(raw) {
  const parsed = parseDominantThresholdInput(raw);
  if (!Number.isFinite(parsed)) {
    document.getElementById('dominantThresholdVal').value = formatDominantThreshold(state.dominantThreshold);
    toast('Invalid threshold', true);
    return;
  }
  state.dominantThreshold = Math.max(DOMINANT_THRESHOLD_MIN, Math.min(DOMINANT_THRESHOLD_MAX, parsed));
  document.getElementById('dominantThresholdSlider').value = dominantThresholdToSlider(state.dominantThreshold).toFixed(3);
  document.getElementById('dominantThresholdVal').value = formatDominantThreshold(state.dominantThreshold);
  renderDominantHistogram();
}

function handleDominantThresholdKey(e) {
  if (e.key === 'Enter') e.currentTarget.blur();
  if (e.key === 'Escape') {
    e.currentTarget.value = formatDominantThreshold(state.dominantThreshold);
    e.currentTarget.blur();
  }
}

function resetDominantZoom() {
  state.dominantZoom.start = 1;
  state.dominantZoom.end = null;
  state.dominantZoom.barWidth = null;
  renderDominantHistogram();
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
      const flippedCount = state.orientationStats.flipped.length;
      document.getElementById('runSummary').innerHTML =
        `Positions scored: <strong style="color:var(--accent)">${L}</strong><br>` +
        `Candidates found: <strong style="color:var(--accent)">${kept.length}</strong><br>` +
        `Reverse-complement fixed: <strong style="color:var(--accent)">${flippedCount}</strong><br>` +
        `Score range: <strong style="color:var(--accent)">${kept.length
          ? kept[kept.length - 1].score.toFixed(3) + ' – ' + kept[0].score.toFixed(3) : '—'
        }</strong><br>` +
        `Window size: <strong style="color:var(--accent)">±${state.WINDOW} nt</strong>`;
      renderDominantHistogram();

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
  if (!document.getElementById('seqScroll')) return;

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
    return `<div class="candidate-item" id="cand-${idx}" onclick="selectCandidate(${idx})" onmouseenter="highlightCandidateBars(${idx})" onmouseleave="clearCandidateBarHighlight()">
      <span class="cand-pos">${c.start + 1}–${c.start + c.len}</span>
      <span class="cand-len">${c.len}nt</span>
      <span class="cand-seq">${getPrimerSeq(c)}</span>
      <div class="score-bar-wrap"><div class="score-bar" style="width:${(c.score * 100).toFixed(0)}%;background:${sc};"></div></div>
      <span class="cand-score" style="color:${sc};">${c.score.toFixed(3)}</span>
    </div>`;
  }).join('');
}

function setHistogramCandidateHighlight(cand) {
  document.querySelectorAll('.dominant-bar.primer-hover').forEach(el => el.classList.remove('primer-hover'));
  if (!cand) return;
  const lo = cand.start + 1;
  const hi = cand.start + cand.len;
  document.querySelectorAll('.dominant-bar').forEach(el => {
    const pos = parseInt(el.dataset.pos, 10);
    if (pos >= lo && pos <= hi) el.classList.add('primer-hover');
  });
}

function highlightCandidateBars(idx) {
  state.hoveredCandidate = state.candidates[idx] || null;
  setHistogramCandidateHighlight(state.hoveredCandidate);
}

function clearCandidateBarHighlight() {
  state.hoveredCandidate = null;
  setHistogramCandidateHighlight(null);
}

function baseColor(base) {
  return { A: '#147a43', T: '#a8333c', G: '#a58900', C: '#2e609f' }[base] || 'var(--text)';
}

function probToLogAxis(prob) {
  const p = Math.max(0, Math.min(0.999999, prob));
  return -Math.log10(1 - p);
}

function logAxisToProb(value) {
  return Math.max(0, Math.min(0.999999, 1 - Math.pow(10, -value)));
}

function dominantHeight(prob, threshold, useLog) {
  if (prob < threshold) return 0;
  const lo = probToLogAxis(threshold);
  const hi = probToLogAxis(DOMINANT_THRESHOLD_MAX);
  const scaled = (probToLogAxis(prob) - lo) / Math.max(0.001, hi - lo);
  return Math.max(0, Math.min(1, scaled)) * 100;
}

function formatYAxisProbability(value, isTop = false) {
  if (isTop) return '1';
  if (value >= 0.9995) return value.toFixed(7).replace(/0+$/, '').replace(/\.$/, '');
  return value.toFixed(3).replace(/0+$/, '').replace(/\.$/, '');
}

function updateDominantYAxis(threshold) {
  const loGap = Math.max(DOMINANT_THRESHOLD_MAX_GAP, 1 - threshold);
  const hiGap = DOMINANT_THRESHOLD_MAX_GAP;
  const tickProb = frac => 1 - loGap * Math.pow(hiGap / loGap, frac);
  const marks = [
    ['dominantY100', 1, true],
    ['dominantY80', tickProb(0.8), false],
    ['dominantY60', tickProb(0.6), false],
    ['dominantY40', tickProb(0.4), false],
    ['dominantY20', tickProb(0.2), false],
    ['dominantY0', threshold, false],
  ];
  marks.forEach(([id, value, isTop]) => {
    const el = document.getElementById(id);
    if (el) el.textContent = formatYAxisProbability(value, isTop);
  });
}

function dominantBreakdown(posIndex) {
  const col = state.colDistribs[posIndex];
  if (!col) return [];
  return ['A', 'T', 'G', 'C'].map(base => ({
    base,
    count: col.cnt[base] || 0,
    prob: col.freq[base] || 0,
  })).sort((a, b) => b.prob - a.prob);
}

function renderDominantHistogram() {
  const histEl = document.getElementById('dominantHist');
  const posEl = document.getElementById('dominantPositions');
  const xAxisEl = document.getElementById('dominantXAxis');
  const tip = document.getElementById('dominantTooltip');
  const scroller = document.getElementById('dominantChartScroll');
  const panel = document.getElementById('uploadZone');
  if (!histEl || !posEl) return;
  if (!state.dominantStats.length) {
    panel?.classList.add('no-data');
    histEl.style.removeProperty('--bar-width');
    histEl.innerHTML = `<div class="dominant-empty">
      <div class="upload-icon">⊕</div>
      <div class="upload-title">Drop FASTA here</div>
      <div class="upload-sub">Dominant base probabilities render immediately after upload</div>
    </div>`;
    if (xAxisEl) xAxisEl.innerHTML = '';
    posEl.textContent = 'No analysis yet';
    return;
  }
  panel?.classList.remove('no-data');

  const threshold = state.dominantThreshold;
  updateDominantYAxis(threshold);
  const useLog = false;
  const zoomStart = Math.max(1, state.dominantZoom.start || 1);
  const zoomEnd = Math.min(state.dominantStats.length, state.dominantZoom.end || state.dominantStats.length);
  const visibleStats = state.dominantStats.slice(zoomStart - 1, zoomEnd);
  const above = visibleStats.filter(d => d.prob >= threshold);
  const barWidth = resolveDominantBarWidth(visibleStats.length);
  const barGap = barWidth <= 3 ? 0 : 3;
  histEl.style.setProperty('--bar-width', `${barWidth}px`);
  histEl.style.setProperty('--bar-gap', `${barGap}px`);
  if (xAxisEl) {
    xAxisEl.style.setProperty('--bar-width', `${barWidth}px`);
    xAxisEl.style.setProperty('--bar-gap', `${barGap}px`);
  }
  histEl.innerHTML = visibleStats.map(d => {
    const isAbove = d.prob >= threshold;
    const inHover = state.hoveredCandidate && d.pos >= state.hoveredCandidate.start + 1 && d.pos <= state.hoveredCandidate.start + state.hoveredCandidate.len;
    const h = dominantHeight(d.prob, threshold, useLog);
    return `<span class="dominant-bar ${isAbove ? 'above' : 'below'} ${inHover ? 'primer-hover' : ''}" data-pos="${d.pos}"
      style="height:${h.toFixed(1)}%;background:${baseColor(d.base)};"
      title="Pos ${d.pos}: ${d.base} ${(d.prob*100).toFixed(2)}% (${d.count}/${d.nvalid})"><span>${d.base}</span></span>`;
  }).join('');
  if (xAxisEl) {
    xAxisEl.innerHTML = visibleStats.map(d =>
      `<span class="dominant-x-tick">${d.pos === 1 || d.pos % 5 === 0 ? d.pos : ''}</span>`
    ).join('');
  }

  histEl.querySelectorAll('.dominant-bar').forEach(bar => {
    const pos = parseInt(bar.dataset.pos, 10);
    const stat = state.dominantStats[pos - 1];
    bar.addEventListener('mouseenter', () => {
      if (!tip || !stat) return;
      const rows = dominantBreakdown(pos - 1).map(v =>
        `<div class="dominant-tip-row"><span class="dominant-tip-swatch" style="background:${baseColor(v.base)}"></span><span>${v.base} ${(v.prob*100).toFixed(1)}%</span></div>`
      ).join('');
      tip.innerHTML = `<b>Pos ${pos}</b> <span style="color:#aaa;">n=${stat.nvalid}</span>${rows}`;
      tip.style.display = 'block';
    });
    bar.addEventListener('mousemove', e => {
      if (!tip) return;
      tip.style.left = `${Math.min(e.clientX + 14, window.innerWidth - tip.offsetWidth - 12)}px`;
      tip.style.top = `${Math.max(e.clientY - tip.offsetHeight - 12, 12)}px`;
    });
    bar.addEventListener('mouseleave', () => {
      if (tip) tip.style.display = 'none';
    });
  });

  posEl.textContent = `View ${zoomStart}-${zoomEnd} · ${above.length}/${visibleStats.length} positions >= ${threshold.toFixed(3)}`;
  if (scroller) scroller.scrollLeft = Math.min(scroller.scrollLeft, scroller.scrollWidth);
}

function resolveDominantBarWidth(nBars) {
  if (state.dominantZoom.barWidth) return state.dominantZoom.barWidth;
  const scroller = document.getElementById('dominantChartScroll');
  const available = Math.max(1, scroller?.clientWidth || 1);
  if (!nBars) return 1;
  return Math.max(1, Math.floor(available / nBars));
}

function getDominantBarGap(nBars) {
  return resolveDominantBarWidth(nBars) <= 3 ? 0 : 3;
}

function saveDominantPng() {
  if (!state.dominantStats.length) { toast('Run analysis before saving chart', true); return; }

  const useLog = false;
  const threshold = state.dominantThreshold;
  const zoomStart = Math.max(1, state.dominantZoom.start || 1);
  const zoomEnd = Math.min(state.dominantStats.length, state.dominantZoom.end || state.dominantStats.length);
  const visibleStats = state.dominantStats.slice(zoomStart - 1, zoomEnd);
  const barW = Math.max(8, Math.min(32, state.dominantZoom.barWidth || 18));
  const gap = 3;
  const margin = { left: 70, right: 22, top: 36, bottom: 62 };
  const plotH = 220;
  const width = Math.max(960, margin.left + margin.right + visibleStats.length * (barW + gap));
  const height = margin.top + plotH + margin.bottom;
  const canvas = document.createElement('canvas');
  canvas.width = width;
  canvas.height = height;
  const ctx = canvas.getContext('2d');

  ctx.fillStyle = '#0d0f14';
  ctx.fillRect(0, 0, width, height);
  ctx.strokeStyle = '#2a3045';
  ctx.fillStyle = '#7a8aaa';
  ctx.font = '12px Arial';
  ctx.textAlign = 'right';
  for (let pct = 0; pct <= 100; pct += 20) {
    const y = margin.top + plotH - (pct / 100) * plotH;
    const loGap = Math.max(DOMINANT_THRESHOLD_MAX_GAP, 1 - threshold);
    const hiGap = DOMINANT_THRESHOLD_MAX_GAP;
    const labelValue = pct === 100 ? 1 : 1 - loGap * Math.pow(hiGap / loGap, pct / 100);
    ctx.beginPath();
    ctx.moveTo(margin.left, y);
    ctx.lineTo(width - margin.right, y);
    ctx.stroke();
    ctx.fillText(formatYAxisProbability(labelValue, pct === 100), margin.left - 8, y + 4);
  }

  visibleStats.forEach((d, idx) => {
    const isAbove = d.prob >= threshold;
    const h = dominantHeight(d.prob, threshold, useLog) / 100 * plotH;
    const x = margin.left + idx * (barW + gap);
    const y = margin.top + plotH - h;
    ctx.fillStyle = isAbove ? baseColor(d.base) : '#3b4250';
    ctx.fillRect(x, y, barW, h);
    if (isAbove && h > 18) {
      ctx.fillStyle = '#d0dcea';
      ctx.font = 'bold 13px monospace';
      ctx.textAlign = 'center';
      ctx.fillText(d.base, x + barW / 2, y + 14);
    }
    if (d.pos === 1 || d.pos % 5 === 0) {
      ctx.fillStyle = '#7a8aaa';
      ctx.font = '12px Arial';
      ctx.fillText(d.pos, x + barW / 2, margin.top + plotH + 18);
    }
  });

  ctx.fillStyle = '#d0dcea';
  ctx.font = '14px Arial';
  ctx.textAlign = 'left';
  ctx.fillText(`Positions ${zoomStart}-${zoomEnd} · Threshold: ${threshold.toFixed(3)}${useLog ? ' · log scale' : ''}`, margin.left, 22);
  ctx.textAlign = 'right';
  ctx.fillText('Position', width - margin.right, height - 22);

  const link = document.createElement('a');
  link.download = 'dominant-base-probability.png';
  link.href = canvas.toDataURL('image/png');
  link.click();
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
    cell.addEventListener('mouseenter', () => {
      const varLines = vd.sorted.slice(0, 5).map(v =>
        `  <span style="color:${baseColor(v.base)};">${v.base}</span>` +
        `  ${v.count.toString().padStart(5)}  <b style="color:var(--accent3)">${(v.pct*100).toFixed(2)}%</b>`
      ).join('<br>');
      const gapLine = vd.ngap > 0
        ? `<br>  <span style="color:var(--text-faint);">-</span>  ${vd.ngap.toString().padStart(5)}  <span style="color:var(--text-faint)">${(vd.ngap/total*100).toFixed(2)}%</span>`
        : '';
      tip.innerHTML =
        `<b style="color:var(--text)">Pos ${ls.pos}</b>  base: <b style="color:${baseColor(base)};">${base}</b><br>` +
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
      for (let i = 0; i < Math.max(qSeq.length, sSeq.length); i++) {
        const q = qSeq[i] || '-', s = sSeq[i] || '-';
        const same = q === s, gap = q === '-' || s === '-';
        identLine += same ? '|' : (gap ? ' ' : '.');
        function bspan(b, mis) {
          if (b === '-') return '<span style="color:var(--text-faint)">-</span>';
          const bc = baseColor(b);
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
      const c = baseColor(b);
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
    const parsedSeqs = parseFASTA(e.target.result);
    if (!parsedSeqs.length) { toast('No valid FASTA sequences found', true); return; }

    const { seqs, flipped } = orientSequencesToReference(parsedSeqs);
    const lens = [...new Set(seqs.map(s => s.seq.length))];
    let seqsAligned = seqs;
    if (lens.length > 1) {
      toast(`Warning: unequal lengths. Using first ${lens[0]} nt.`);
      seqsAligned = seqs.map(s => ({ ...s, seq: s.seq.slice(0, lens[0]) }));
    }

    state.sequences      = seqsAligned;
    state.orientationStats = { flipped };
    state.seqLen         = seqsAligned[0].seq.length;
    state.consensusSeq   = buildConsensus(seqsAligned);
    state.candidates     = [];
    const { raw, windowed, colDistribs } = computeJSD(seqsAligned, state.WINDOW, state.LAMBDA_JSD);
    state.jsdScores      = raw;
    state.windowedScores = windowed;
    state.colDistribs    = colDistribs;
    state.dominantStats  = dominantBaseStats(colDistribs);
    state.activeCandidate = null;
    state.hoveredCandidate = null;
    state.showFull       = false;

    document.getElementById('uploadZone').style.display    = 'flex';
    document.getElementById('nSeqDisplay').textContent     = seqsAligned.length;
    document.getElementById('lenDisplay').textContent      = `${state.seqLen} nt`;
    document.getElementById('seqLabel').className          = 'panel-label active';
    document.getElementById('runBtn').disabled             = false;
    document.getElementById('candidatesList').innerHTML = '<div class="no-candidates">Run analysis to find candidates</div>';
    document.getElementById('candidateCount').textContent = '';
    document.getElementById('statsSection').style.display      = 'none';

    resetDominantZoom();
    setStatus(`Loaded ${seqsAligned.length} sequences × ${state.seqLen} nt · histogram ready${flipped.length ? ` · ${flipped.length} RC fixed` : ''}`, 'blue');
    toast(`Loaded: ${seqsAligned.length} sequences${flipped.length ? ` · fixed ${flipped.length} reverse-complement` : ''}`);
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

function bindDominantSelection() {
  const scroller = document.getElementById('dominantChartScroll');
  const selection = document.getElementById('dominantSelection');
  if (!scroller || !selection) return;

  let dragStartX = null;
  let dragStartScroll = 0;

  function visibleStats() {
    const zoomStart = Math.max(1, state.dominantZoom.start || 1);
    const zoomEnd = Math.min(state.dominantStats.length, state.dominantZoom.end || state.dominantStats.length);
    return state.dominantStats.slice(zoomStart - 1, zoomEnd);
  }

  function pointerToIndex(clientX) {
    const rect = scroller.getBoundingClientRect();
    const x = clientX - rect.left + scroller.scrollLeft;
    const slot = resolveDominantBarWidth(visibleStats().length) + getDominantBarGap(visibleStats().length);
    const idx = Math.floor(Math.max(0, x) / slot);
    return Math.max(0, Math.min(visibleStats().length - 1, idx));
  }

  scroller.addEventListener('pointerdown', e => {
    if (!state.dominantStats.length || e.button !== 0) return;
    if (e.detail > 1) return;
    dragStartX = e.clientX;
    dragStartScroll = scroller.scrollLeft;
    const rect = scroller.getBoundingClientRect();
    selection.style.display = 'block';
    selection.style.left = `${e.clientX - rect.left + scroller.scrollLeft}px`;
    selection.style.width = '1px';
    scroller.setPointerCapture(e.pointerId);
  });

  scroller.addEventListener('pointermove', e => {
    if (dragStartX === null) return;
    const rect = scroller.getBoundingClientRect();
    const current = e.clientX - rect.left + scroller.scrollLeft;
    const start = dragStartX - rect.left + dragStartScroll;
    selection.style.left = `${Math.max(0, Math.min(start, current))}px`;
    selection.style.width = `${Math.abs(current - start)}px`;
  });

  scroller.addEventListener('pointerup', e => {
    if (dragStartX === null) return;
    const visible = visibleStats();
    const moved = Math.abs(e.clientX - dragStartX);
    const startIdx = pointerToIndex(dragStartX);
    const endIdx = pointerToIndex(e.clientX);
    dragStartX = null;
    selection.style.display = 'none';

    if (moved < 8 || !visible.length) return;
    const lo = Math.min(startIdx, endIdx);
    const hi = Math.max(startIdx, endIdx);
    if (hi <= lo) return;

    state.dominantZoom.start = visible[lo].pos;
    state.dominantZoom.end = visible[hi].pos;
    renderDominantHistogram();
    scroller.scrollLeft = 0;
  });

  scroller.addEventListener('pointercancel', () => {
    dragStartX = null;
    selection.style.display = 'none';
  });

  scroller.addEventListener('dblclick', e => {
    if (!state.dominantStats.length) return;
    e.preventDefault();
    resetDominantZoom();
    scroller.scrollLeft = 0;
  });

  scroller.addEventListener('wheel', e => {
    if (!state.dominantStats.length) return;
    e.preventDefault();
    const oldWidth = state.dominantZoom.barWidth || resolveDominantBarWidth(visibleStats().length);
    const nextWidth = Math.max(1, Math.min(80, oldWidth + (e.deltaY < 0 ? 2 : -2)));
    if (nextWidth === oldWidth) return;

    const rect = scroller.getBoundingClientRect();
    const anchor = e.clientX - rect.left + scroller.scrollLeft;
    const ratio = anchor / Math.max(1, scroller.scrollWidth);
    state.dominantZoom.barWidth = nextWidth;
    renderDominantHistogram();
    scroller.scrollLeft = Math.max(0, ratio * scroller.scrollWidth - (e.clientX - rect.left));
  }, { passive: false });
}

function bindLayoutResizers() {
  const root = document.documentElement;
  const main = document.querySelector('.main-layout');
  const mainSplitter = document.getElementById('mainSplitter');
  const candidateSplitter = document.getElementById('candidateSplitter');
  const leftPanel = document.getElementById('leftPanel');
  const candidates = document.getElementById('candidatesSection');
  if (!main || !mainSplitter || !candidateSplitter || !leftPanel || !candidates) return;

  mainSplitter.addEventListener('pointerdown', e => {
    e.preventDefault();
    mainSplitter.classList.add('dragging');
    mainSplitter.setPointerCapture(e.pointerId);

    function move(ev) {
      const rect = main.getBoundingClientRect();
      const rightWidth = Math.max(280, Math.min(rect.width - 340, rect.right - ev.clientX));
      root.style.setProperty('--right-width', `${rightWidth}px`);
      renderDominantHistogram();
    }
    function up() {
      mainSplitter.classList.remove('dragging');
      mainSplitter.removeEventListener('pointermove', move);
      mainSplitter.removeEventListener('pointerup', up);
      mainSplitter.removeEventListener('pointercancel', up);
    }
    mainSplitter.addEventListener('pointermove', move);
    mainSplitter.addEventListener('pointerup', up);
    mainSplitter.addEventListener('pointercancel', up);
  });

  candidateSplitter.addEventListener('pointerdown', e => {
    e.preventDefault();
    candidateSplitter.classList.add('dragging');
    candidateSplitter.setPointerCapture(e.pointerId);

    function move(ev) {
      const rect = leftPanel.getBoundingClientRect();
      const height = Math.max(120, Math.min(rect.height - 180, rect.bottom - ev.clientY));
      root.style.setProperty('--candidates-height', `${height}px`);
      renderDominantHistogram();
    }
    function up() {
      candidateSplitter.classList.remove('dragging');
      candidateSplitter.removeEventListener('pointermove', move);
      candidateSplitter.removeEventListener('pointerup', up);
      candidateSplitter.removeEventListener('pointercancel', up);
    }
    candidateSplitter.addEventListener('pointermove', move);
    candidateSplitter.addEventListener('pointerup', up);
    candidateSplitter.addEventListener('pointercancel', up);
  });
}

function bindChartResizeObserver() {
  const scroller = document.getElementById('dominantChartScroll');
  if (!scroller || typeof ResizeObserver === 'undefined') return;
  let lastWidth = scroller.clientWidth;
  const observer = new ResizeObserver(() => {
    if (Math.abs(scroller.clientWidth - lastWidth) < 2) return;
    lastWidth = scroller.clientWidth;
    if (!state.dominantZoom.barWidth) renderDominantHistogram();
  });
  observer.observe(scroller);
}

// Init sliders from CONFIG
syncK(document.getElementById('kSlider').value);
syncC(document.getElementById('cSlider').value);
document.getElementById('dominantThresholdSlider').value = dominantThresholdToSlider(state.dominantThreshold).toFixed(3);
syncDominantThreshold(document.getElementById('dominantThresholdSlider').value);
bindDominantSelection();
bindLayoutResizers();
bindChartResizeObserver();
