import { useState, useCallback } from "react";

/* ═══ CODON TABLE ═══ */
const CT={TTT:"F",TTC:"F",TTA:"L",TTG:"L",CTT:"L",CTC:"L",CTA:"L",CTG:"L",ATT:"I",ATC:"I",ATA:"I",ATG:"M",GTT:"V",GTC:"V",GTA:"V",GTG:"V",TCT:"S",TCC:"S",TCA:"S",TCG:"S",CCT:"P",CCC:"P",CCA:"P",CCG:"P",ACT:"T",ACC:"T",ACA:"T",ACG:"T",GCT:"A",GCC:"A",GCA:"A",GCG:"A",TAT:"Y",TAC:"Y",TAA:"*",TAG:"*",CAT:"H",CAC:"H",CAA:"Q",CAG:"Q",GAT:"D",GAC:"D",GAA:"E",GAG:"E",TGA:"*",TGT:"C",TGC:"C",TGG:"W",CGT:"R",CGC:"R",CGA:"R",CGG:"R",AGT:"S",AGC:"S",AGA:"R",AGG:"R",GGT:"G",GGC:"G",GGA:"G",GGG:"G"};
const toAA=c=>CT[c]||"?";
const rc=s=>s.split("").reverse().map(b=>({A:"T",T:"A",G:"C",C:"G"}[b]||"N")).join("");

/* ═══ GENBANK PARSER (handles \r\n) ═══ */
function parseGB(rawText) {
  const text = rawText.replace(/\r\n/g, "\n").replace(/\r/g, "\n");
  const feats = []; let seq = ""; let inS = false; let inF = false; let cur = null; let lq = "";
  for (const ln of text.split("\n")) {
    if (ln.startsWith("ORIGIN")) { inS = true; inF = false; continue; }
    if (ln.startsWith("//")) break;
    if (ln.startsWith("FEATURES")) { inF = true; continue; }
    if (inS) { seq += ln.replace(/[^a-zA-Z]/g, ""); continue; }
    if (inF) {
      if (ln.length > 5 && ln.charAt(5) !== " " && ln.charAt(0) === " ") {
        if (cur) feats.push(cur);
        const parts = ln.trim().split(/\s+/);
        cur = { type: parts[0], loc: parts.slice(1).join(""), q: {} }; lq = "";
      } else if (cur && ln.match(/^\s{21}\//)) {
        const qm = ln.trim().match(/^\/(\w+)=?"?([^"]*)"?$/);
        if (qm) { cur.q[qm[1]] = qm[2]; lq = qm[1]; }
      } else if (cur && ln.match(/^\s{21}[^/]/) && !cur.loc.includes(")")) {
        cur.loc += ln.trim();
      } else if (cur && lq && ln.match(/^\s{21}/)) {
        cur.q[lq] = (cur.q[lq] || "") + ln.trim().replace(/"/g, "");
      }
    }
  }
  if (cur) feats.push(cur);
  return { seq: seq.toUpperCase(), feats };
}

function getCDS(gb) {
  for (const f of gb.feats) {
    if (f.type !== "CDS") continue;
    const nums = f.loc.match(/\d+/g);
    if (!nums || nums.length < 2) continue;
    const segs = [];
    for (let i = 0; i < nums.length - 1; i += 2)
      segs.push([parseInt(nums[i]) - 1, parseInt(nums[i + 1])]);
    let cds = "";
    for (const [s, e] of segs) cds += gb.seq.slice(s, e);
    const name = gb.feats.find(x => x.type === "gene")?.q?.label || "Gene";
    return { segs, cds, prot: Math.floor(cds.length / 3) - 1, name };
  }
  return null;
}

/* ═══ COORDINATE HELPERS ═══ */
function c2g(segs, cp) { let c = 0; for (const [s, e] of segs) { if (c + (e - s) > cp) return s + (cp - c); c += (e - s); } return null; }
function g2aa(segs, gp) { let c = 0; for (const [s, e] of segs) { if (s <= gp && gp < e) return Math.floor((c + gp - s) / 3) + 1; c += (e - s); } return null; }
function gCodon(segs, seq, n) { const g = c2g(segs, (n - 1) * 3); if (g === null) return null; return { g, cod: seq.slice(g, g + 3), aa: toAA(seq.slice(g, g + 3)) }; }

/* ═══ gRNA FINDER ═══ */
function findG(seq, tgt, rng = 50) {
  const r = [];
  for (let p = Math.max(0, tgt - rng - 23); p <= Math.min(seq.length - 23, tgt + rng); p++) {
    const pam = seq.slice(p + 20, p + 23);
    if (["AGG", "TGG", "CGG", "GGG"].includes(pam)) {
      const sp = seq.slice(p, p + 20), cut = p + 17;
      const gc = Math.round(sp.split("").filter(b => b === "G" || b === "C").length / 20 * 100);
      if (!sp.includes("TTTTT") && !sp.includes("AAAAA") && Math.abs(cut - tgt) <= rng)
        r.push({ sp, pam, str: "+", cut, d: cut - tgt, gc, ps: p });
    }
    const pr = seq.slice(p, p + 3);
    if (["CCA", "CCT", "CCC", "CCG"].includes(pr)) {
      const sp = rc(seq.slice(p + 3, p + 23)), pam2 = rc(pr), cut = p + 6;
      const gc = Math.round(sp.split("").filter(b => b === "G" || b === "C").length / 20 * 100);
      if (!sp.includes("TTTTT") && !sp.includes("AAAAA") && Math.abs(cut - tgt) <= rng)
        r.push({ sp, pam: pam2, str: "-", cut, d: cut - tgt, gc, ps: p });
    }
  }
  return r;
}

function armType(g, mp) { return g.str === "+" ? (mp > g.cut ? "PROX" : "DIST") : (mp < g.cut ? "PROX" : "DIST"); }

/* ═══ SILENT MUTATION FINDER ═══ */
function findSilent(seq, segs, g) {
  let pamS;
  if (g.str === "+") pamS = g.ps + 20;
  else { const t = rc(g.sp); const i = seq.indexOf(t); if (i < 0) return null; pamS = i - 3; }
  if (pamS < 0 || pamS + 3 > seq.length) return null;
  const pis = g.str === "+" ? [1, 2] : [0, 1];
  for (const pi of pis) {
    const gp = pamS + pi; if (gp < 0 || gp >= seq.length) continue;
    const an = g2aa(segs, gp); if (!an) continue;
    const cg = c2g(segs, (an - 1) * 3); if (cg === null) continue;
    const cod = seq.slice(cg, cg + 3), cp = gp - cg, oa = toAA(cod);
    for (const alt of ["A", "C", "G", "T"]) {
      if (alt === cod[cp]) continue;
      const mc = cod.slice(0, cp) + alt + cod.slice(cp + 1);
      if (toAA(mc) === oa) {
        const tp = seq.slice(pamS, pamS + 3).split(""); tp[pi] = alt;
        const ok = g.str === "+" ? tp.slice(1).join("") !== "GG" : rc(tp.join("")).slice(1) !== "GG";
        if (ok) { const po = g.str === "+" ? seq.slice(pamS, pamS + 3) : rc(seq.slice(pamS, pamS + 3)); const pn = g.str === "+" ? tp.join("") : rc(tp.join(""));
          return { gp, nb: alt, lb: `p.${oa}${an}${oa}`, oc: cod, nc: mc, pur: `PAM ${po}→${pn}` }; }
      }
    }
  }
  for (let si = 10; si < 20; si++) {
    let gp; if (g.str === "+") gp = g.ps + si;
    else { const t = rc(g.sp); const i = seq.indexOf(t); if (i < 0) continue; gp = i + 19 - si; }
    if (gp === undefined || gp < 0 || gp >= seq.length) continue;
    const an = g2aa(segs, gp); if (!an) continue;
    const cg = c2g(segs, (an - 1) * 3); if (cg === null) continue;
    const cod = seq.slice(cg, cg + 3), cp = gp - cg, oa = toAA(cod);
    for (const alt of ["A", "C", "G", "T"]) { if (alt === cod[cp]) continue; const mc = cod.slice(0, cp) + alt + cod.slice(cp + 1);
      if (toAA(mc) === oa) return { gp, nb: alt, lb: `p.${oa}${an}${oa}`, oc: cod, nc: mc, pur: `Seed pos ${si + 1}/20` }; }
  }
  return null;
}

/* ═══ ssODN BUILDER ═══ */
function mkODN(seq, g, mps, mbs, sils = []) {
  let dS, dE;
  if (g.str === "+") { dS = g.cut - 36; dE = g.cut + 91; } else { dS = g.cut - 91; dE = g.cut + 36; }
  if (dS < 0 || dE > seq.length) return null;
  const pl = seq.slice(dS, dE).split(""), wt = seq.slice(dS, dE);
  mps.forEach((m, i) => { const idx = m - dS; if (idx >= 0 && idx < 127) pl[idx] = mbs[i]; });
  sils.forEach(s => { const idx = s.gp - dS; if (idx >= 0 && idx < 127) pl[idx] = s.nb; });
  const od = g.str === "+" ? rc(pl.join("")) : pl.join("");
  const wo = g.str === "+" ? rc(wt) : wt;
  const df = []; for (let i = 0; i < 127; i++) if (od[i] !== wo[i]) df.push(i);
  return { od, wo, df, sl: g.str === "+" ? "−strand (target)" : "+strand (target)" };
}

/* ═══ TAG INSERT SEQUENCES ═══ */
const TAGS = {
  "SD40-2xHA": { ins: "GCTAAAGCCAAAAACAACCAGGGATCCGGACTGTTGCTGTTCTGCCCTATTTGCGGGTTTACATGTCGCCAGAAGGGCAACTTACTTCGCCATATTAACCTGCACACAGGGGAAAAGTTATTTAAGTACCACCTGTATGGCGGCTACCCCTACGACGTGCCCGACTACGCCGGCTATCCGTATGATGTCCCGGACTATGCATAA", d: "Linker(30)+SD40(108)+2×HA(60)+STOP = 204bp" },
  "2xHA": { ins: "GCTAAAGCCAAAAACAACCAGGGATCCGGAGGCGGCTACCCCTACGACGTGCCCGACTACGCCGGCTATCCGTATGATGTCCCGGACTATGCATAA", d: "Linker(30)+2×HA(60)+STOP = 96bp" },
};

/* ═══ POINT MUTATION DESIGN ═══ */
function designPM(gb, mutStr) {
  const cds = getCDS(gb); if (!cds) return { err: "No CDS annotation found in GenBank file" };
  const S = gb.seq;
  const m = mutStr.match(/([A-Z])(\d+)([A-Z])/i);
  if (!m) return { err: "Cannot parse mutation. Use format like L72S, E280D, R176C" };
  const [, wA, ns, mA] = m; const an = parseInt(ns);
  const ci = gCodon(cds.segs, S, an);
  if (!ci) return { err: `Cannot map amino acid ${an} to genomic position. Check numbering.` };
  if (ci.aa !== wA.toUpperCase()) return { err: `Expected ${wA.toUpperCase()} at position ${an} but found ${ci.aa} (${ci.cod}). Verify the amino acid numbering matches the transcript.` };
  let bM = null, bC = [];
  for (const b1 of "ACGT") for (const b2 of "ACGT") for (const b3 of "ACGT") {
    const mc = b1 + b2 + b3; if (toAA(mc) === mA.toUpperCase()) {
      const ch = []; for (let i = 0; i < 3; i++) if (ci.cod[i] !== mc[i]) ch.push({ p: ci.g + i, w: ci.cod[i], m: mc[i] });
      if (!bM || ch.length < bC.length) { bM = mc; bC = ch; }
    }
  }
  if (!bM) return { err: `No codon encodes ${mA.toUpperCase()}` };
  const mps = bC.map(c => c.p), mbs = bC.map(c => c.m);
  const ag = findG(S, ci.g, 50), pg = ag.filter(g => armType(g, ci.g) === "PROX").sort((a, b) => Math.abs(a.d) - Math.abs(b.d));
  if (!pg.length) return { err: "No gRNAs found that place mutation in the 91nt proximal arm. The region may be too AT-rich for SpCas9." };
  const sel = pg.slice(0, 2);
  const res = { type: "pm", gene: cds.name, an, wA: wA.toUpperCase(), mA: mA.toUpperCase(), wC: ci.cod, mC: bM, gp: ci.g, ch: bC, gs: [], os: [], ss: [], ps: [] };
  sel.forEach((g, i) => { try {
    const sm = findSilent(S, cds.segs, g), sa = sm ? [sm] : []; const od = mkODN(S, g, mps, mbs, sa);
    res.gs.push({ n: `gRNA${i + 1}`, sp: g.sp, pm: g.pam, str: g.str, gc: g.gc, d: g.d, arm: `91nt proximal ✓ (${Math.abs(g.d)}bp)` });
    if (od) res.os.push({ ...od, n: `ssODN${i + 1}`, gi: i }); if (sm) res.ss.push({ ...sm, gi: i + 1 });
  } catch(e) { console.error('gRNA design error:', e); }});
  const fwP = Math.max(0, ci.g - 250);
  res.ps = [{ n: "Fw", s: S.slice(fwP, fwP + 24) }, { n: "Rev", s: rc(S.slice(Math.min(S.length - 24, ci.g + 250), Math.min(S.length, ci.g + 274))) }];
  return res;
}

/* ═══ C-TERMINAL KI DESIGN ═══ */
function designCT(gb, tag, haL) {
  const cds = getCDS(gb); if (!cds) return { err: "No CDS annotation found in GenBank file" };
  const S = gb.seq;
  const ls = cds.segs[cds.segs.length - 1];
  const ss = ls[1] - 3;
  const sc = S.slice(ss, ss + 3);
  if (!["TAA", "TAG", "TGA"].includes(sc)) return { err: `Expected stop codon at end of CDS, found "${sc}". Check CDS annotation.` };
  const ti = TAGS[tag]; if (!ti) return { err: `Tag "${tag}" not available. Use SD40-2xHA or 2xHA.` };
  const hl = parseInt(haL) || 250;
  const h5s = Math.max(0, ss - hl), h3e = Math.min(S.length, ss + 3 + hl);
  const h5 = S.slice(h5s, ss), h3 = S.slice(ss + 3, h3e);
  const ag = findG(S, ss, 50).sort((a, b) => Math.abs(a.d) - Math.abs(b.d));
  if (!ag.length) return { err: "No SpCas9 gRNAs found within 50bp of stop codon. Region may be too AT-rich." };
  const sel = ag.slice(0, 2);
  const sils = [];
  sel.forEach((g, i) => { try {
    const sm = findSilent(S, cds.segs, g);
    if (sm && sm.gp >= h5s && sm.gp < ss) sils.push({ ...sm, gi: i + 1 });
  } catch(e) { /* silent mutation search failed, proceed without */ }});
  const h5a = h5.split("");
  sils.forEach(s => { const idx = s.gp - h5s; if (idx >= 0 && idx < h5a.length) h5a[idx] = s.nb; });
  const donor = h5a.join("") + ti.ins + h3;
  const la = gCodon(cds.segs, S, cds.prot);
  return { type: "ct", gene: cds.name, stop: sc, sp: ss + 1, prot: cds.prot,
    lastAA: la ? `${la.aa}${cds.prot}` : "?", tag, td: ti.d, il: ti.ins.length,
    hl, h5l: h5.length, h3l: h3.length, dl: donor.length, donor,
    gs: sel.map((g, i) => ({ n: `gRNA${i + 1}`, sp: g.sp, pm: g.pam, str: g.str, gc: g.gc, d: g.d,
      note: `Cut ${Math.abs(g.d)}bp ${g.d < 0 ? "5′" : "3′"} of stop` })),
    ss: sils,
    ps: [{ n: "Fw", s: S.slice(Math.max(0, h5s - 50), Math.max(0, h5s - 50) + 24) },
         { n: "Rev", s: rc(S.slice(Math.min(S.length - 24, h3e + 25), Math.min(S.length, h3e + 49))) }],
    amp: `WT ~${h3e + 49 - h5s + 50}bp · KI ~${h3e + 49 - h5s + 50 + ti.ins.length}bp` };
}

/* ═══ KNOCKOUT DESIGN ═══ */
function designKO(gb) {
  const cds = getCDS(gb); if (!cds) return { err: "No CDS annotation found" };
  const S = gb.seq;
  const tSegs = cds.segs.slice(1, 4);
  if (!tSegs.length) return { err: "Not enough coding exons for KO design" };
  let best = tSegs[0], bl = 0;
  tSegs.forEach(([s, e]) => { if (e - s > bl) { bl = e - s; best = [s, e]; } });
  const [exS, exE] = best, exM = Math.floor((exS + exE) / 2);
  const si = cds.segs.findIndex(s => s[0] === exS) + 1;
  const ag = findG(S, exM, Math.floor(bl / 2) + 10).filter(g => g.cut >= exS + 5 && g.cut <= exE - 5);
  if (!ag.length) return { err: "No gRNAs found within target exon" };
  ag.sort((a, b) => b.gc - a.gc);
  return { type: "ko", gene: cds.name, exon: `CDS segment ${si} (${exS + 1}–${exE}, ${bl}bp)`,
    exSz: bl, prot: cds.prot,
    gs: ag.slice(0, 2).map((g, i) => ({ n: `gRNA${i + 1}`, sp: g.sp, pm: g.pam, str: g.str, gc: g.gc, d: g.d,
      note: `Cut at ${g.cut + 1}, ${g.cut - exS}bp into exon` })),
    ps: [{ n: "Fw", s: S.slice(Math.max(0, exS - 200), Math.max(0, exS - 200) + 24) },
         { n: "Rev", s: rc(S.slice(Math.min(S.length - 24, exE + 175), Math.min(S.length, exE + 199))) }],
    amp: `~${exE + 199 - exS + 200}bp`,
    strat: "NHEJ-mediated frameshift via Cas9 RNP. Two gRNAs for redundancy. Indels cause premature stop → NMD. Screen by Sanger + ICE/TIDE. Confirm protein loss by Western blot." };
}

/* ═══════════════════════ UI ═══════════════════════ */
const K = { bg: "#070d18", s1: "#0e1628", s2: "#162038", bd: "#1e2f4d", ac: "#22d3ee", a2: "#818cf8", wn: "#fbbf24", er: "#fb7185", ok: "#34d399", t: "#e2e8f0", tm: "#8b9bb5", td: "#4b5c73" };

export default function App() {
  const [pt, sPt] = useState("");
  const [gbRaw, sGb] = useState("");
  const [fn, sFn] = useState("");
  const [mut, sMut] = useState("");
  const [tag, sTag] = useState("SD40-2xHA");
  const [ha, sHa] = useState("250");
  const [res, sRes] = useState(null);
  const [err, sErr] = useState("");
  const [dbg, sDbg] = useState("");

  const onF = useCallback(e => {
    const f = e.target.files[0]; if (!f) return;
    sFn(f.name); sRes(null); sErr("");
    const r = new FileReader();
    r.onload = ev => { sGb(ev.target.result); sDbg(`File loaded: ${f.name} (${ev.target.result.length} chars)`); };
    r.onerror = () => sErr("Failed to read file");
    r.readAsText(f);
  }, []);

  const run = () => {
    sErr(""); sRes(null); sDbg("");
    if (!gbRaw) { sErr("Upload a GenBank file first"); return; }
    if (!pt) { sErr("Select a project type first"); return; }
    try {
      const gb = parseGB(gbRaw);
      sDbg(d => d + `\nParsed: ${gb.seq.length}bp, ${gb.feats.length} features`);
      if (!gb.seq) { sErr("Could not parse DNA sequence from file. Is this a GenBank (.gb) file?"); return; }
      const cds = getCDS(gb);
      if (!cds) { sErr("No CDS annotation found. The GenBank file needs a CDS feature with join(...) coordinates."); return; }
      sDbg(d => d + `\nCDS: ${cds.name}, ${cds.segs.length} segments, ${cds.prot} aa`);
      let r;
      if (pt === "pm") { if (!mut) { sErr("Enter a mutation (e.g. L72S)"); return; } r = designPM(gb, mut); }
      else if (pt === "ct") { r = designCT(gb, tag, ha); }
      else { r = designKO(gb); }
      if (r.err) { sErr(r.err); return; }
      sRes(r);
    } catch (e) { sErr(`Error: ${e.message}`); console.error(e); }
  };

  const Bd = ({ children, color = K.ac }) => (
    <span style={{ fontSize: 9, fontWeight: 700, padding: "2px 6px", borderRadius: 3, background: `${color}15`, color, border: `1px solid ${color}33`, whiteSpace: "nowrap" }}>{children}</span>
  );

  return (
    <div style={{ minHeight: "100vh", background: K.bg, color: K.t, fontFamily: "system-ui,-apple-system,sans-serif" }}>
      <div style={{ borderBottom: `1px solid ${K.bd}`, padding: "10px 16px", display: "flex", alignItems: "center", gap: 8, flexWrap: "wrap" }}>
        <div style={{ width: 26, height: 26, borderRadius: 5, background: `linear-gradient(135deg,${K.ac},${K.a2})`, display: "flex", alignItems: "center", justifyContent: "center", fontSize: 12, color: "#000", fontWeight: 800 }}>A</div>
        <b style={{ fontSize: 13 }}>ASSURED CRISPR Designer</b>
        <Bd>⚡ Instant</Bd>
      </div>

      <div style={{ maxWidth: 880, margin: "0 auto", padding: "16px 12px" }}>
        {/* Project type */}
        <div style={{ display: "flex", gap: 6, marginBottom: 12, flexWrap: "wrap" }}>
          {[{ id: "pm", ic: "🎯", lb: "Point Mutation", sb: "ssODN · 127nt · ASSURED" },
            { id: "ct", ic: "🏷️", lb: "C-terminal Tag KI", sb: "SD40 · 2xHA · donor" },
            { id: "ko", ic: "✂️", lb: "Gene Knockout", sb: "Frameshift · NHEJ" }
          ].map(p => (
            <div key={p.id} onClick={() => sPt(p.id)} style={{ flex: 1, minWidth: 140, padding: "10px 12px", borderRadius: 8, cursor: "pointer",
              border: `1.5px solid ${pt === p.id ? K.ac : K.bd}`, background: pt === p.id ? `${K.ac}08` : K.s1 }}>
              <div style={{ fontSize: 16 }}>{p.ic}</div>
              <div style={{ fontWeight: 700, fontSize: 11, color: pt === p.id ? K.ac : K.t }}>{p.lb}</div>
              <div style={{ fontSize: 9, color: K.td }}>{p.sb}</div>
            </div>
          ))}
        </div>

        {/* Inputs */}
        <div style={{ display: "flex", gap: 6, marginBottom: 10, flexWrap: "wrap", alignItems: "center" }}>
          <label style={{ flex: "1 1 180px", display: "flex", alignItems: "center", gap: 6, padding: "9px 10px", borderRadius: 7,
            border: `1.5px dashed ${fn ? K.ok : K.bd}`, background: K.s1, cursor: "pointer" }}>
            <input type="file" accept=".gb,.gbk,.genbank,.txt" onChange={onF} style={{ display: "none" }} />
            <span style={{ fontSize: 14 }}>{fn ? "✅" : "📁"}</span>
            <span style={{ fontSize: 11, color: fn ? K.ok : K.tm }}>{fn || "Upload .gb file"}</span>
          </label>
          {pt === "pm" && <input value={mut} onChange={e => sMut(e.target.value)} placeholder="e.g. L72S"
            onKeyDown={e => { if (e.key === "Enter") run(); }}
            style={{ width: 110, padding: "9px 10px", borderRadius: 7, border: `1px solid ${K.bd}`, background: K.s1, color: K.t, fontSize: 13, fontWeight: 600 }} />}
          {pt === "ct" && <>
            <select value={tag} onChange={e => sTag(e.target.value)}
              style={{ padding: 9, borderRadius: 7, border: `1px solid ${K.bd}`, background: K.s1, color: K.t, fontSize: 11 }}>
              {Object.keys(TAGS).map(t => <option key={t} value={t}>{t}</option>)}
            </select>
            <select value={ha} onChange={e => sHa(e.target.value)}
              style={{ padding: 9, borderRadius: 7, border: `1px solid ${K.bd}`, background: K.s1, color: K.t, fontSize: 11 }}>
              <option value="250">250bp HA</option><option value="500">500bp HA</option><option value="750">750bp HA</option>
            </select>
          </>}
          <button onClick={run} style={{ padding: "9px 20px", borderRadius: 7, border: "none", fontWeight: 700, fontSize: 12,
            background: `linear-gradient(135deg,${K.ac},${K.a2})`, color: "#000", cursor: "pointer" }}>⚡ Design</button>
        </div>

        {err && <div style={{ padding: 8, borderRadius: 7, background: `${K.er}12`, border: `1px solid ${K.er}33`, color: K.er, fontSize: 11, marginBottom: 10 }}>{err}</div>}

        {/* ═══ RESULTS ═══ */}
        {res && (
          <div style={{ background: K.s1, borderRadius: 10, border: `1px solid ${K.bd}`, overflow: "hidden" }}>
            <div style={{ padding: "12px 14px", borderBottom: `1px solid ${K.bd}`, display: "flex", alignItems: "center", gap: 8, flexWrap: "wrap" }}>
              <span style={{ fontSize: 18 }}>{res.type === "pm" ? "🎯" : res.type === "ct" ? "🏷️" : "✂️"}</span>
              <div>
                <div style={{ fontWeight: 700, fontSize: 14 }}>
                  {res.type === "pm" && `${res.gene} p.${res.wA}${res.an}${res.mA}`}
                  {res.type === "ct" && `${res.gene} C-terminal ${res.tag}`}
                  {res.type === "ko" && `${res.gene} Knockout`}
                </div>
                <div style={{ fontSize: 10, color: K.tm }}>
                  {res.type === "pm" && `${res.wC}→${res.mC} · ${res.ch.length} base change${res.ch.length > 1 ? "s" : ""} · Pos ${res.gp + 1}`}
                  {res.type === "ct" && `Stop ${res.stop} at ${res.sp} · Last aa: ${res.lastAA} · ${res.tag} (${res.td}) · Donor: ${res.dl}bp`}
                  {res.type === "ko" && `Target: ${res.exon} · ${res.prot} aa protein`}
                </div>
              </div>
              <div style={{ marginLeft: "auto" }}><Bd color={K.ok}>{res.type === "pm" ? "ASSURED ✓" : res.type === "ct" ? "KI Design ✓" : "KO Design ✓"}</Bd></div>
            </div>

            <div style={{ padding: 14 }}>
              {/* gRNAs */}
              <div style={{ fontSize: 10, fontWeight: 700, color: K.ac, marginBottom: 6, textTransform: "uppercase", letterSpacing: 1 }}>gRNA Sequences</div>
              {res.gs.map((g, i) => (
                <div key={i} style={{ display: "flex", alignItems: "center", gap: 5, padding: "6px 8px", background: K.s2, borderRadius: 5, marginBottom: 3, border: `1px solid ${K.bd}`, flexWrap: "wrap", fontSize: 11 }}>
                  <b style={{ minWidth: 44 }}>{g.n}</b>
                  <code style={{ fontSize: 11, color: K.ac, letterSpacing: 0.5 }}>{g.sp}</code>
                  <b>{g.pm}</b>
                  <Bd color={g.str === "+" ? K.ac : K.a2}>{g.str}str</Bd>
                  <span style={{ color: K.tm, fontSize: 9 }}>GC={g.gc}%</span>
                  {g.arm && <Bd color={K.ok}>{g.arm}</Bd>}
                  {g.note && <span style={{ color: K.td, fontSize: 9 }}>{g.note}</span>}
                </div>
              ))}

              {/* ssODNs (PM) */}
              {res.os?.map((s, i) => (
                <div key={i} style={{ marginTop: 12 }}>
                  <div style={{ fontSize: 9, fontWeight: 700, color: K.a2, marginBottom: 4, textTransform: "uppercase", letterSpacing: 1 }}>{s.n} · {s.sl} · 127nt (36+91)</div>
                  <div style={{ background: K.s2, borderRadius: 6, padding: 8, border: `1px solid ${K.bd}`, overflowX: "auto" }}>
                    {[{ l: "WT   ", q: s.wo, w: true }, { l: "ssODN", q: s.od, w: false }].map(({ l, q, w }) => (
                      <div key={l} style={{ whiteSpace: "nowrap", marginBottom: w ? 1 : 0, lineHeight: 1.4 }}>
                        <span style={{ fontFamily: "monospace", fontSize: 9, color: K.td, fontWeight: 700 }}>{l} 5′-</span>
                        {q.split("").map((b, j) => (
                          <span key={j} style={{ fontFamily: "monospace", fontSize: 10,
                            color: s.df.includes(j) ? K.er : (w ? K.td : K.t), fontWeight: s.df.includes(j) ? 800 : 400,
                            textDecoration: w && s.df.includes(j) ? "line-through" : "none",
                            background: !w && s.df.includes(j) ? `${K.wn}22` : "transparent" }}>{b}</span>
                        ))}
                        <span style={{ fontFamily: "monospace", fontSize: 9, color: K.td }}>-3′</span>
                      </div>
                    ))}
                    <div style={{ marginTop: 4, fontSize: 9, color: K.ok, fontWeight: 600 }}>✓ 36nt distal + 91nt proximal = 127nt · Mutation in proximal arm · {s.df.length} changes</div>
                  </div>
                </div>
              ))}

              {/* Donor (CT KI) */}
              {res.type === "ct" && res.donor && (
                <div style={{ marginTop: 12 }}>
                  <div style={{ fontSize: 9, fontWeight: 700, color: K.a2, marginBottom: 4, textTransform: "uppercase", letterSpacing: 1 }}>
                    Donor · {res.dl}bp ({res.h5l} 5′HA + {res.il} insert + {res.h3l} 3′HA)
                  </div>
                  <div style={{ background: K.s2, borderRadius: 6, padding: 8, border: `1px solid ${K.bd}`, overflowX: "auto", maxHeight: 180, overflowY: "auto" }}>
                    {Array.from({ length: Math.ceil(res.donor.length / 80) }, (_, li) => {
                      const st = li * 80;
                      return (
                        <div key={li} style={{ whiteSpace: "nowrap", lineHeight: 1.3 }}>
                          <span style={{ fontFamily: "monospace", fontSize: 9, color: K.td }}>{String(st + 1).padStart(5)} </span>
                          {res.donor.slice(st, st + 80).split("").map((b, j) => {
                            const p = st + j; let cl;
                            if (p < res.h5l) cl = "#3b82f6";
                            else if (p < res.h5l + 30) cl = "#8b5cf6";
                            else if (p < res.h5l + res.il - 63) cl = "#10b981";
                            else if (p < res.h5l + res.il - 3) cl = "#f59e0b";
                            else if (p < res.h5l + res.il) cl = "#ef4444";
                            else cl = "#0891b2";
                            return <span key={j} style={{ fontFamily: "monospace", fontSize: 9, color: cl }}>{b}</span>;
                          })}
                        </div>
                      );
                    })}
                  </div>
                  <div style={{ display: "flex", gap: 6, marginTop: 4, flexWrap: "wrap", fontSize: 8 }}>
                    {[["5′HA", "#3b82f6"], ["Linker", "#8b5cf6"], ["Tag", "#10b981"], ["2×HA", "#f59e0b"], ["STOP", "#ef4444"], ["3′HA", "#0891b2"]].map(([l, c]) => (
                      <span key={l} style={{ color: c, fontWeight: 700 }}>■ {l}</span>
                    ))}
                  </div>
                </div>
              )}

              {/* KO strategy */}
              {res.type === "ko" && (
                <div style={{ marginTop: 12, padding: 8, background: K.s2, borderRadius: 6, border: `1px solid ${K.bd}` }}>
                  <div style={{ fontSize: 9, fontWeight: 700, color: K.wn, marginBottom: 3, textTransform: "uppercase" }}>Strategy</div>
                  <div style={{ fontSize: 11, color: K.tm, lineHeight: 1.5 }}>{res.strat}</div>
                </div>
              )}

              {/* Silent mutations */}
              {res.ss?.length > 0 && (
                <div style={{ marginTop: 12 }}>
                  <div style={{ fontSize: 9, fontWeight: 700, color: K.wn, marginBottom: 4, textTransform: "uppercase", letterSpacing: 1 }}>Silent Mutations</div>
                  {res.ss.map((m, i) => (
                    <div key={i} style={{ fontSize: 10, padding: "3px 6px", background: `${K.wn}08`, borderRadius: 3, border: `1px solid ${K.wn}22`, marginBottom: 2 }}>
                      <b style={{ color: K.wn }}>gRNA{m.gi}:</b> <span style={{ color: K.tm }}>{m.lb} ({m.oc}→{m.nc}) — {m.pur}</span>
                    </div>
                  ))}
                </div>
              )}

              {/* Primers */}
              <div style={{ marginTop: 12 }}>
                <div style={{ fontSize: 9, fontWeight: 700, color: K.ac, marginBottom: 4, textTransform: "uppercase", letterSpacing: 1 }}>Primers{res.amp ? ` · ${res.amp}` : ""}</div>
                {res.ps.map((p, i) => (
                  <div key={i} style={{ display: "flex", gap: 6, padding: "4px 6px", background: K.s2, borderRadius: 3, marginBottom: 2, border: `1px solid ${K.bd}`, fontSize: 10 }}>
                    <b style={{ minWidth: 24 }}>{p.n}</b>
                    <code style={{ color: K.tm, fontSize: 10 }}>{p.s}</code>
                  </div>
                ))}
              </div>
            </div>
          </div>
        )}
      </div>
    </div>
  );
}
