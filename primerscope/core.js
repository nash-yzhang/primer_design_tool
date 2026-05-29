// ── CONFIG ──────────────────────────────────────────────────────────────────
// Default parameters (mirrors config.json)
const CONFIG = {
  k:          0.60,
  c:          0.70,
  dominantThreshold: 0.900,
  minLen:     18,
  maxLen:     24,
  window:     3,
  lambdaJSD:  0.5,
  trimLimit:  500,
};

const BASES = ['A', 'T', 'G', 'C'];

// ── FASTA PARSING ────────────────────────────────────────────────────────────
function parseFASTA(text) {
  const seqs = [];
  let cur = null;
  for (const line of text.split('\n')) {
    const l = line.trim();
    if (!l) continue;
    if (l.startsWith('>')) {
      if (cur) seqs.push(cur);
      cur = { name: l.slice(1).split(/\s+/)[0], seq: '' };
    } else if (cur) {
      cur.seq += l.toUpperCase().replace(/[^ATGCN\-\.]/g, '-');
    }
  }
  if (cur) seqs.push(cur);
  return seqs;
}

function reverseComplement(seq) {
  const comp = { A: 'T', T: 'A', G: 'C', C: 'G', N: 'N', '-': '-', '.': '.' };
  let out = '';
  for (let i = seq.length - 1; i >= 0; i--) out += comp[seq[i]] || 'N';
  return out;
}

function alignedIdentity(refSeq, querySeq) {
  const L = Math.min(refSeq.length, querySeq.length);
  let matches = 0;
  let compared = 0;
  for (let i = 0; i < L; i++) {
    const r = refSeq[i];
    const q = querySeq[i];
    if (!BASES.includes(r) || !BASES.includes(q)) continue;
    compared++;
    if (r === q) matches++;
  }
  return { identity: compared ? matches / compared : 0, compared };
}

function orientSequencesToReference(seqs, minDelta = 0.05) {
  if (seqs.length < 2) return { seqs, flipped: [] };

  const refSeq = seqs[0].seq;
  const flipped = [];
  const oriented = seqs.map((s, idx) => {
    if (idx === 0) return { ...s, orientation: 'reference' };

    const forward = alignedIdentity(refSeq, s.seq);
    const rcSeq = reverseComplement(s.seq);
    const reverse = alignedIdentity(refSeq, rcSeq);
    const shouldFlip = reverse.compared >= 8 && reverse.identity > forward.identity + minDelta;

    if (!shouldFlip) return { ...s, orientation: 'forward' };

    flipped.push({
      name: s.name,
      forwardIdentity: forward.identity,
      reverseIdentity: reverse.identity,
    });
    return { ...s, seq: rcSeq, orientation: 'reverse-complement' };
  });

  return { seqs: oriented, flipped };
}

function buildConsensus(seqs) {
  if (!seqs.length) return '';
  const L = seqs[0].seq.length;
  let consensus = '';
  for (let i = 0; i < L; i++) {
    const cnt = {};
    for (const s of seqs) {
      const b = s.seq[i] || '-';
      if (b !== '-' && b !== '.') cnt[b] = (cnt[b] || 0) + 1;
    }
    const entries = Object.entries(cnt).sort((a, b) => b[1] - a[1]);
    consensus += entries.length ? entries[0][0] : '-';
  }
  return consensus;
}

// ── JSD COMPUTATION ──────────────────────────────────────────────────────────
function entropy(p) {
  return -p.reduce((s, pi) => (pi <= 0 ? s : s + pi * Math.log2(pi)), 0);
}

function jsdColumn(colFreqs, bgFreq, lambda) {
  const bases = ['A', 'T', 'G', 'C'];
  const r   = bases.map(b => lambda * (colFreqs[b] || 0) + (1 - lambda) * (bgFreq[b] || 0.25));
  const Hr  = entropy(r);
  const Hpc = entropy(bases.map(b => colFreqs[b] || 0));
  const Hq  = entropy(bases.map(b => bgFreq[b]   || 0.25));
  return Math.max(0, Math.min(1, Hr - (lambda * Hpc + (1 - lambda) * Hq)));
}

function computeJSD(seqs, windowSize, lambda) {
  const L = seqs[0].seq.length;
  const bgFreq = { A: 0.25, T: 0.25, G: 0.25, C: 0.25 };
  const raw = [];
  const colDistribs = [];

  for (let i = 0; i < L; i++) {
    const cnt = { A: 0, T: 0, G: 0, C: 0 };
    let ngap = 0;
    for (const s of seqs) {
      const b = s.seq[i] || '-';
      if (b === '-' || b === '.') ngap++;
      else if (Object.prototype.hasOwnProperty.call(cnt, b)) cnt[b]++;
    }
    const nvalid  = seqs.length - ngap;
    const gapFrac = ngap / seqs.length;

    const freq = { A: 0, T: 0, G: 0, C: 0 };
    if (nvalid > 0) {
      for (const b of BASES) freq[b] = cnt[b] / nvalid;
    }
    colDistribs.push({ freq, cnt: { ...cnt }, nvalid, gapFrac });

    let score = 0;
    if (gapFrac <= 0.7) {
      score = jsdColumn(freq, bgFreq, lambda) * (1 - gapFrac);
    }
    raw.push(score);
  }

  // Window smoothing: 0.5 * raw(i) + 0.5 * mean(window around i)
  const W = windowSize;
  const windowed = [];
  for (let i = 0; i < L; i++) {
    const lo = Math.max(0, i - W);
    const hi = Math.min(L - 1, i + W);
    let sum = 0, n = 0;
    for (let j = lo; j <= hi; j++) { sum += raw[j]; n++; }
    windowed.push(0.5 * raw[i] + 0.5 * (sum / n));
  }

  // Normalize to [0, 1]
  const mx = Math.max(...windowed, 1e-9);
  const normed = windowed.map(v => v / mx);

  return { raw, windowed: normed, colDistribs };
}

function dominantBaseStats(colDistribs) {
  return colDistribs.map((col, idx) => {
    const entries = BASES.map(base => ({
      base,
      prob: col.freq[base] || 0,
      count: col.cnt[base] || 0,
    })).sort((a, b) => b.prob - a.prob);
    return { pos: idx + 1, ...entries[0], nvalid: col.nvalid, gapFrac: col.gapFrac };
  });
}

// 3'-weighted primer suitability score
function primerScore(winScores, start, len) {
  let weightedSum = 0, totalWeight = 0;
  for (let i = 0; i < len; i++) {
    const w = i + 1; // increases toward 3' end
    weightedSum += w * winScores[start + i];
    totalWeight += w;
  }
  return weightedSum / totalWeight;
}
