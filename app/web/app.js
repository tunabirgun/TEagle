"use strict";
const $ = s => document.querySelector(s);
const css = n => getComputedStyle(document.documentElement).getPropertyValue(n).trim();
const state = { seq: "", pair: null, structural: [], elemSpan: null, source: null, pcrPairs: [] };

/* ---------- source citations (clickable, verified DOIs — mirror backend refs.py) ---------- */
const REFLINKS = {
  Wicker2007:   { url:"https://doi.org/10.1038/nrg2165", cite:"Wicker T, et al. (2007) A unified classification system for eukaryotic transposable elements. Nat Rev Genet 8:973-982." },
  Pfam:         { url:"https://www.ebi.ac.uk/interpro/", cite:"Mistry J, et al. (2021) Pfam: the protein families database in 2021. Nucleic Acids Res 49:D412-D419." },
  HMMER:        { url:"https://doi.org/10.1371/journal.pcbi.1002195", cite:"Eddy SR (2011) Accelerated Profile HMM Searches. PLoS Comput Biol 7:e1002195." },
  Dfam:         { url:"https://doi.org/10.1186/s13100-020-00230-y", cite:"Storer J, et al. (2021) The Dfam community resource of transposable element families. Mob DNA 12:2." },
  RepeatMasker: { url:"https://www.repeatmasker.org/", cite:"Smit AFA, Hubley R, Green P. RepeatMasker Open-4.0." },
  NCBI:         { url:"https://www.ncbi.nlm.nih.gov/nuccore/", cite:"NCBI Entrez / E-utilities (Sayers E, NCBI)." },
  Primer3:      { url:"https://doi.org/10.1093/nar/gks596", cite:"Untergasser A, et al. (2012) Primer3 — new capabilities and interfaces. Nucleic Acids Res 40:e115." },
  minimap2:     { url:"https://doi.org/10.1093/bioinformatics/bty191", cite:"Li H (2018) Minimap2: pairwise alignment for nucleotide sequences. Bioinformatics 34:3094-3100." },
};
function srcChip(key, url){
  const r = REFLINKS[key]; if(!r) return "";
  return `<a class="srcchip" href="${url||r.url}" target="_blank" rel="noopener" title="Source — ${r.cite}" onclick="event.stopPropagation()">source ↗</a>`;
}
/* plain-language glossary for table headers (hover to learn the abbreviation) */
const GLOSSARY = {
  "Feature":"Structural hallmark found in the sequence — e.g. LTR, TIR, target-site duplication, poly-A tail.",
  "Segment":"An annotated gene part from the record: exon (kept in mRNA), intron (spliced out), or CDS (protein-coding).",
  "Splice site":"The two bases at each intron boundary (donor…acceptor); canonical eukaryotic introns are GT…AG (or GC–AG / AT–AC).",
  "Intron span (0-based)":"Intron location in the loaded sequence, 0-based half-open [start, end).",
  "Span (0-based)":"Location in the sequence, 0-based half-open [start, end).",
  "Coords (0-based)":"Location in the sequence, 0-based half-open [start, end).",
  "Coords":"Location of this amplicon in the searched sequence.",
  "Len":"Length in base pairs.",
  "Metric":"Feature-specific measure — terminal-repeat identity %, a motif, or a length.",
  "Method":"How TEagle detected this feature (the algorithm/heuristic used).",
  "aa":"Length of the predicted protein (open reading frame) in amino acids.",
  "ORF":"Open reading frame — a stretch translatable with no internal stop codon (≥40 aa here).",
  "Class/family":"TE class and superfamily (Wicker 2007 scheme), e.g. LTR/Copia.",
  "Dfam family":"The specific named family in the Dfam library, e.g. Copia_I or L1HS.",
  "Str":"Strand of the match — + forward, − reverse complement.",
  "Div":"Divergence — % difference between your sequence and the Dfam family consensus (lower = closer).",
  "Score":"RepeatMasker alignment score — higher means stronger homology.",
  "Pair":"Which designed primer pair produced this amplicon.",
  "Source":"The sequence that was searched — your specimen or a custom background.",
  "Mism F/R":"Mismatches in the forward / reverse primer binding site (the 3′ end is kept exact).",
  "Call":"On-target = amplicon at the intended locus; off-target = amplified elsewhere.",
};
function glossify(root){
  if(!root) return;
  root.querySelectorAll("th").forEach(th => {
    const key = th.textContent.trim();
    if(GLOSSARY[key] && !th.title){ th.title = GLOSSARY[key]; th.classList.add("gloss"); }
  });
}
/* calibrated ETA for heavy jobs (annotate / splice) — accuracy improves as this PC logs runs */
function fmtEta(s){ s = Math.round(s); return s >= 90 ? `~${Math.round(s/60)} min` : `~${Math.max(1,s)}s`; }
async function etaText(job, size){
  try {
    const e = await api("/api/eta", { job, size });
    const basis = e.basis === "calibrated" ? ` (from ${e.n} past run${e.n>1?"s":""} on this PC)` : " (first run — rough estimate)";
    return `ETA ${fmtEta(e.eta_s)}${basis}`;
  } catch { return ""; }
}
function liveTimer(host){
  const t0 = Date.now();
  const h = setInterval(() => { const c = host.querySelector(".etaclock"); if(c) c.textContent = `${((Date.now()-t0)/1000).toFixed(0)}s elapsed`; }, 500);
  return () => clearInterval(h);
}

/* ---------- deterministic real sample element (seeded, so structure detection fires) ---------- */
function mulberry32(a){return function(){a|=0;a=a+0x6D2B79F5|0;let t=Math.imul(a^a>>>15,1|a);t=t+Math.imul(t^t>>>7,61|t)^t;return ((t^t>>>14)>>>0)/4294967296;};}
const B="ACGT", STOP={TAA:1,TAG:1,TGA:1};
function mkSample(){
  const r=mulberry32(20260717);
  const rb=n=>{let s="";for(let i=0;i<n;i++)s+=B[(r()*4)|0];return s;};
  const rc=()=>{let c;do{c=B[(r()*4)|0]+B[(r()*4)|0]+B[(r()*4)|0];}while(STOP[c]);return c;};
  const orf=nc=>{let s="ATG";for(let i=0;i<nc;i++)s+=rc();return s+"TAA";};
  const mut=(s,rate)=>{const a=s.split("");for(let i=0;i<a.length;i++)if(r()<rate)a[i]=B[(r()*4)|0];return a.join("");};
  const tsd=rb(5), ltr=rb(160), internal=rb(300)+orf(210)+rb(280);
  const seq=rb(250)+tsd+ltr+internal+mut(ltr,0.02)+tsd+rb(250);
  return ">sample_LTR_element  demo construct (illustrative)\n"+seq.replace(/(.{60})/g,"$1\n");
}

/* ---------- api ---------- */
async function api(path, body){
  const r = await fetch(path, {method:"POST", headers:{"Content-Type":"application/json"}, body:JSON.stringify(body)});
  const j = await r.json();
  if(!r.ok) throw new Error(j.error || ("HTTP "+r.status));
  return j;
}
function err(msg){ $("#errbox").innerHTML = `<div class="err-banner">⚠ ${msg}</div>`; }
function clearErr(){ $("#errbox").innerHTML = ""; }

/* ---------- click-to-FASTA submenu (any TE feature) + helpers ---------- */
const _RC = {A:"T",T:"A",G:"C",C:"G",U:"A",R:"Y",Y:"R",S:"S",W:"W",K:"M",M:"K",B:"V",D:"H",H:"D",V:"B",N:"N"};
function revcomp(s){ return s.toUpperCase().split("").reverse().map(c=>_RC[c]||"N").join(""); }
function cleanSeq(){ return (state.seq || (seqEl&&seqEl.value) || "").replace(/^>.*$/gm,"").replace(/[^A-Za-z]/g,"").toUpperCase(); }
function slice(s,e){ return cleanSeq().slice(s,e); }
let _menu, _menuTrigger;
function closeMenu(restore){ if(_menu){ _menu.remove(); _menu=null; const t=_menuTrigger; _menuTrigger=null; if(restore && t && t.focus) t.focus(); } }
document.addEventListener("click", e => { if(_menu && !_menu.contains(e.target)) closeMenu(); });
document.addEventListener("keydown", e => { if(e.key==="Escape" && _menu) closeMenu(true); });
/* keyboard activation for clickable rows / domain rows (Enter / Space) */
document.addEventListener("keydown", e => {
  if((e.key === "Enter" || e.key === " ") && document.activeElement){
    const el = document.activeElement;
    if(el.classList && (el.classList.contains("clickable") || el.classList.contains("domrow"))){ e.preventDefault(); el.click(); }
  }
});
const _reduce = () => matchMedia("(prefers-reduced-motion: reduce)").matches;
function scrollToCard(el){ if(el) el.scrollIntoView({ behavior: _reduce() ? "auto" : "smooth", block: "start" }); }
function expandCard(el, scroll){ if(!el) return; el.classList.remove("collapsed"); const ch=el.querySelector(".ch"); if(ch) ch.setAttribute("aria-expanded","true"); if(scroll) scrollToCard(el); }

/* ================= publication-quality SVG figures (export PNG/SVG, customizable bg) ================= */
const OK = { RT:"#0072B2", INT:"#E69F00", RNaseH:"#009E73", PR:"#CC79A7", GAG:"#7A7A7A", CHR:"#D55E00", TPase:"#D55E00",
             LTR:"#56B4E9", TIR:"#E69F00", tail:"#CC79A7", ORF:"#4C6C97", on:"#009E73", off:"#D55E00", ladder:"#999999" };
const _figState = {}; // per-figure background choice
function _paper(bg){ return bg==="white" ? "#FFFFFF" : bg==="transparent" ? "none" : "#0C1116"; }
function _ink(bg){ return bg==="dark" ? "#E6EDF1" : "#333333"; }        // transparent -> dark ink (reads on white/manuscript)
function _faint(bg){ return bg==="dark" ? "#8A959D" : "#6b7075"; }
function _grid(bg){ return bg==="dark" ? "#26313B" : "#cfd3d7"; }
function labelInk(hex){                                                 // pick black/white label by best WCAG contrast
  const lin = v => { v/=255; return v<=0.03928 ? v/12.92 : Math.pow((v+0.055)/1.055, 2.4); };
  const L = 0.2126*lin(parseInt(hex.slice(1,3),16)) + 0.7152*lin(parseInt(hex.slice(3,5),16)) + 0.0722*lin(parseInt(hex.slice(5,7),16));
  return (1.05/(L+0.05)) >= ((L+0.05)/0.05) ? "#fff" : "#111";
}
const esc = s => String(s).replace(/&/g,"&amp;").replace(/</g,"&lt;").replace(/>/g,"&gt;");
const FIGFONT = "Arial, Helvetica, sans-serif";   // publication-convention figure font
const GENECOL = { exon:"#009E73", intron:"#8792a0", cds:"#D55E00" };   // exon/CDS Okabe-Ito; intron neutral connector

/* ================= genome-browser figure viewer (GDV-style: ruler, tracks, wheel-zoom, pan) ================= */
function gvNiceStep(span, ticks){
  const raw = Math.max(1, span/ticks), mag = Math.pow(10, Math.floor(Math.log10(raw))), n = raw/mag;
  return (n < 1.5 ? 1 : n < 3.5 ? 2 : n < 7.5 ? 5 : 10) * mag;
}
function gvTracksFromRec(rec){
  const tracks = [], reps = [];
  (rec.structural||[]).forEach(e => {
    if(e.type.startsWith("LTR") || e.type.startsWith("TIR")){
      const col = e.type.startsWith("LTR") ? OK.LTR : OK.TIR;
      [e.five_prime, e.three_prime].forEach(p => { if(p) reps.push({ start:p[0], end:p[1], color:col, label:e.type.split(" ")[0], span:e.element_span, tip:`${e.type} ${p[0]}–${p[1]}` }); });
    } else if(e.pos){ reps.push({ start:e.pos[0], end:e.pos[1], color:OK.tail, label:e.type.split(" ")[0], tip:`${e.type} ${e.pos[0]}–${e.pos[1]}` }); }
  });
  if(reps.length) tracks.push({ name:"terminal repeats", height:20, features:reps });
  const doms = (rec.domains||[]).map(d => ({ start:d.nt[0], end:d.nt[1], color:OK[d.domain]||"#888", label:d.domain, tip:`${d.domain} · ${d.label} · nt ${d.nt[0]}–${d.nt[1]} · score ${d.score}` }));
  if(doms.length) tracks.push({ name:"protein domains", height:22, features:doms });
  const orfs = (rec.orfs||[]).map(o => ({ start:o.start, end:o.end, color:OK.ORF, strand:o.strand, tip:`ORF ${o.strand}${o.frame} · ${o.length_aa} aa` }));
  if(orfs.length) tracks.push({ name:"ORFs (± strand)", height:26, features:orfs, stranded:true });
  return { length: rec.composition.length || 1, tracks };
}
function gvTracksFromGene(gm, len){
  const tracks = [];
  const feat = (gm.exons||[]).map(e => ({ start:e.start, end:e.end, color:GENECOL.exon, label:"exon", tip:`exon ${e.start}–${e.end} (${e.end-e.start} bp)` }))
    .concat((gm.introns||[]).map(i => ({ start:i.start, end:i.end, color:GENECOL.intron, intron:true, tip:`intron ${i.start}–${i.end}${i.donor?` · ${i.donor}…${i.acceptor}${i.canonical?" (canonical)":""}`:""}` })));
  if(feat.length) tracks.push({ name:"exons / introns", height:22, features:feat });
  const cds = (gm.cds||[]).map(c => ({ start:c.start, end:c.end, color:GENECOL.cds, label:"CDS", tip:`CDS ${c.start}–${c.end}` }));
  if(cds.length) tracks.push({ name:"CDS (coding)", height:16, features:cds });
  return { length: len || 1, tracks };
}
function gvTheme(theme, forExport){
  if(forExport) return { paper:"none", ink:"#222", faint:"#555", grid:"#dcdfe3", track:"#00000000", lane:"#0000000d", frame:"#c7ccd2", win:"#1f6feb" };
  if(theme === "white") return { paper:"#ffffff", ink:"#141b21", faint:"#5a6570", grid:"#eceef1", track:"#f6f8fa", lane:"#eef1f4", frame:"#dde1e6", win:"#1f6feb" };
  return { paper:"#0b1016", ink:"#e6edf1", faint:"#8a959d", grid:"#182029", track:"#10171e", lane:"#121b23", frame:"#243039", win:"#4aa8ff" };
}
function svgGenome(model, view, W, theme, forExport){
  const L = model.length, ML = 96, MR = 16, MT = 34, ovH = 13, rulerH = 24;
  const plotW = Math.max(120, W - ML - MR);
  const s0 = view.start, s1 = view.end, span = Math.max(1, s1 - s0);
  const bx = bp => ML + (bp - s0) / span * plotW;
  const ox = bp => ML + bp / L * plotW;
  const T = gvTheme(theme, forExport);
  let y = MT + ovH + 10 + rulerH;
  const trackYs = model.tracks.map(t => { const ty = y; y += (t.height||20) + 20; return ty; });
  const H = y + 12;
  let s = `<svg xmlns="http://www.w3.org/2000/svg" viewBox="0 0 ${W} ${H}" width="${W}" height="${H}" font-family="${FIGFONT}">`;
  s += `<defs><linearGradient id="gvwin" x1="0" x2="0" y1="0" y2="1"><stop offset="0" stop-color="${T.win}" stop-opacity="0.28"/><stop offset="1" stop-color="${T.win}" stop-opacity="0.10"/></linearGradient></defs>`;
  if(T.paper !== "none") s += `<rect width="${W}" height="${H}" fill="${T.paper}"/>`;
  // overview: whole element + current window
  const ovY = MT;
  s += `<text x="${ML}" y="${ovY-5}" fill="${T.faint}" font-size="8.5">whole element · ${L.toLocaleString()} bp</text>`;
  s += `<rect x="${ML}" y="${ovY}" width="${plotW}" height="${ovH}" rx="3" fill="${T.lane}" stroke="${T.frame}"/>`;
  model.tracks.forEach(t => t.features.forEach(f => { const a=ox(f.start), w=Math.max(ox(f.end)-ox(f.start),1); s += `<rect x="${a.toFixed(1)}" y="${ovY+2}" width="${w.toFixed(1)}" height="${ovH-4}" rx="1" fill="${f.color}" opacity="0.6"/>`; }));
  s += `<rect x="${ox(s0).toFixed(1)}" y="${(ovY-2).toFixed(1)}" width="${Math.max(ox(s1)-ox(s0),2).toFixed(1)}" height="${ovH+4}" rx="2" fill="url(#gvwin)" stroke="${T.win}" stroke-width="1.2"/>`;
  // ruler
  const ry = MT + ovH + 10 + rulerH - 7, step = gvNiceStep(span, 7), first = Math.ceil(s0/step)*step;
  s += `<line x1="${ML}" y1="${ry}" x2="${ML+plotW}" y2="${ry}" stroke="${T.frame}"/>`;
  for(let bp=first; bp<=s1+1; bp+=step){ const x=bx(bp); if(x<ML-1||x>ML+plotW+1) continue;
    const lab = bp>=1e6?(bp/1e6)+"M":bp>=1000?(bp/1000)+"k":bp;
    s += `<line x1="${x.toFixed(1)}" y1="${ry}" x2="${x.toFixed(1)}" y2="${ry-5}" stroke="${T.faint}"/>`
      +  `<text x="${x.toFixed(1)}" y="${(ry-8).toFixed(1)}" fill="${T.faint}" font-size="9" text-anchor="middle">${lab}</text>`
      +  `<line x1="${x.toFixed(1)}" y1="${(ry+2).toFixed(1)}" x2="${x.toFixed(1)}" y2="${H-10}" stroke="${T.grid}"/>`; }
  // tracks
  model.tracks.forEach((t,ti) => { const ty=trackYs[ti], th=t.height||20;
    s += `<text x="${ML-10}" y="${(ty+th/2+3).toFixed(1)}" fill="${T.faint}" font-size="9.5" text-anchor="end">${esc(t.name)}</text>`;
    s += `<rect x="${ML}" y="${ty}" width="${plotW}" height="${th}" rx="3" fill="${T.track}"/>`;
    if(t.stranded){ s += `<line x1="${ML}" y1="${(ty+th/2).toFixed(1)}" x2="${ML+plotW}" y2="${(ty+th/2).toFixed(1)}" stroke="${T.grid}"/>`; }
    t.features.forEach(f => { let a=Math.max(bx(f.start),ML), b=Math.min(bx(f.end),ML+plotW); if(b<ML-0.5||a>ML+plotW+0.5) return;
      if(f.intron){ const mid=(a+b)/2; s += `<path d="M ${a.toFixed(1)} ${(ty+th/2).toFixed(1)} L ${mid.toFixed(1)} ${(ty+3).toFixed(1)} L ${b.toFixed(1)} ${(ty+th/2).toFixed(1)}" fill="none" stroke="${f.color}" stroke-width="1.3"><title>${esc(f.tip||"")}</title></path>`; return; }
      const w=Math.max(b-a,1.5), yy = t.stranded ? (f.strand==="+"?ty+2.5:ty+th/2+1.5) : ty+2.5, hh = t.stranded ? th/2-4 : th-5;
      s += `<rect class="gvglyph" x="${a.toFixed(1)}" y="${yy.toFixed(1)}" width="${w.toFixed(1)}" height="${Math.max(hh,3).toFixed(1)}" rx="2.5" fill="${f.color}"><title>${esc(f.tip||"")}</title></rect>`;
      if(f.label && w>26) s += `<text x="${(a+4).toFixed(1)}" y="${(yy+Math.max(hh,3)-3).toFixed(1)}" fill="${labelInk(f.color)}" font-size="9" font-weight="700" pointer-events="none">${esc(f.label)}</text>`; });
  });
  return s + "</svg>";
}
const _gvState = {};
function mountGenomeViewer(host, model, opts){
  if(!host) return; opts = opts || {};
  const key = host.id || opts.base || "gv", base = opts.base || "TEagle_figure", L = model.length;
  const theme = _figState[key+"_th"] || "dark";
  const prev = _gvState[key];
  let view = (prev && prev.len === L) ? { start:prev.start, end:prev.end } : { start:0, end:L };
  host.innerHTML = `<div class="figtoolbar">
      <span class="lbl">bg</span><button class="btn sm" data-th="dark">dark</button><button class="btn sm" data-th="white">light</button>
      <span class="gvpos" data-pos></span><span class="grow"></span>
      <button class="btn sm" data-z="out" aria-label="zoom out">−</button>
      <button class="btn sm" data-z="fit" title="show the whole element">fit</button>
      <button class="btn sm" data-z="in" aria-label="zoom in">+</button>
      <button class="btn sm" data-exp="svg" title="download SVG (transparent background)">⭳ SVG</button>
      <button class="btn sm" data-exp="png" title="download high-res PNG (transparent)">⭳ PNG</button>
    </div>
    <div class="gvview gv-${theme}" data-view tabindex="0" role="group" aria-label="genome viewer — scroll to zoom, drag to pan, arrow keys to navigate, Home to fit">
      <div class="gvsvg" data-svg></div>
      <div class="gvcross" data-cross style="display:none"><span class="gvchip" data-chip></span></div>
    </div>
    <div class="fighint">scroll = zoom (cursor-anchored) · drag = pan · double-click = zoom · arrows / ± / Home = navigate · <b>fit</b> = whole element · exports the current view, transparent</div>`;
  const tb = host.querySelector('[data-th="'+theme+'"]'); if(tb) tb.classList.add("primary");
  const viewEl = host.querySelector("[data-view]"), svgEl = host.querySelector("[data-svg]"), posEl = host.querySelector("[data-pos]"),
        crossEl = host.querySelector("[data-cross]"), chipEl = host.querySelector("[data-chip]");
  const ML = 96, MRr = 16;
  // the SVG is drawn at W = max(clientWidth,320) and CSS-stretched to clientWidth; keep the inverse mapping in
  // sync via the display scale (== 1 at normal widths, corrects the crosshair/click below 320px)
  const _W = () => Math.max(viewEl.clientWidth||620, 320);
  const _scale = () => (viewEl.clientWidth||_W()) / _W();
  const plotW = () => Math.max(120, _W() - ML - MRr);                    // SVG units
  const xToBp = cx => { const r=viewEl.getBoundingClientRect(); const svgX=(cx-r.left)/(_scale()||1);
    const frac=Math.max(0,Math.min(1,(svgX-ML)/plotW())); return view.start + frac*(view.end-view.start); };
  let raf=null, target=null, last=0, lastCX=null, drag=null;
  const save = () => { _gvState[key] = { start:view.start, end:view.end, len:L }; };
  const updateCross = () => { if(lastCX==null || drag){ crossEl.style.display="none"; return; }
    const r=viewEl.getBoundingClientRect(), x=lastCX-r.left, sc=_scale(); if(x < ML*sc || x > (ML+plotW())*sc){ crossEl.style.display="none"; return; }
    crossEl.style.display="block"; crossEl.style.left=x+"px"; chipEl.textContent = Math.round(xToBp(lastCX)).toLocaleString()+" bp"; };
  const draw = () => { svgEl.innerHTML = svgGenome(model, view, Math.max(viewEl.clientWidth||620,320), theme, false);
    if(posEl) posEl.textContent = `${Math.round(view.start).toLocaleString()}–${Math.round(view.end).toLocaleString()} bp · ${((view.end-view.start)/1000).toFixed(2)} kb`;
    updateCross(); save(); };
  const clamp = v => { let sp = Math.min(Math.max(v.end-v.start, 20), L); if(sp>=L){ return {start:0,end:L}; } let st = Math.max(0, Math.min(v.start, L-sp)); return {start:st, end:st+sp}; };
  const animateTo = t => { target = clamp(t); last = 0;
    if(_reduce()){ view = {...target}; if(raf){cancelAnimationFrame(raf);raf=null;} draw(); return; }   // respect reduced-motion
    if(raf) return;
    const tick = now => { const dt = Math.min(now-(last||now), 50); last = now; const e = 1 - Math.exp(-dt/80);   // frame-rate-independent easing
      view.start += (target.start-view.start)*e; view.end += (target.end-view.end)*e;
      if(Math.abs(view.start-target.start)<0.6 && Math.abs(view.end-target.end)<0.6){ view = {...target}; draw(); raf=null; return; }
      draw(); raf = requestAnimationFrame(tick); };
    raf = requestAnimationFrame(tick); };
  const zoomAt = (bp, factor) => { const ref = raf ? target : view; const sp=(ref.end-ref.start)*factor;   // compound on the in-flight target span
    const fracPix=(bp-view.start)/(view.end-view.start); animateTo({ start:bp-fracPix*sp, end:bp+(1-fracPix)*sp }); };
  // wheel: zoom proportional to pixel delta (no trackpad slam), cursor-anchored; ctrlKey = pinch-to-zoom
  viewEl.addEventListener("wheel", e => { e.preventDefault(); const unit = e.deltaMode===1 ? 0.05 : e.deltaMode===2 ? 1 : 0.002;
    zoomAt(xToBp(e.clientX), Math.pow(2, e.deltaY*unit*(e.ctrlKey?10:1))); }, {passive:false});
  // pointerdown: click in the overview band jumps the window; elsewhere starts a pan
  viewEl.addEventListener("pointerdown", e => { const r=viewEl.getBoundingClientRect(), yy=e.clientY-r.top;
    if(yy>=30 && yy<=52){ const sc=_scale(); const frac=Math.max(0,Math.min(1,(e.clientX-r.left-ML*sc)/(plotW()*sc))), bp=frac*L, sp=view.end-view.start; animateTo({start:bp-sp/2, end:bp+sp/2}); return; }
    drag={x:e.clientX, s:view.start, en:view.end}; viewEl.setPointerCapture(e.pointerId); viewEl.classList.add("grabbing"); crossEl.style.display="none"; });
  viewEl.addEventListener("pointermove", e => { lastCX=e.clientX;
    if(drag){ if(raf){cancelAnimationFrame(raf); raf=null; target=null;} const dbp=(e.clientX-drag.x)/plotW()*(drag.en-drag.s); view=clamp({start:drag.s-dbp, end:drag.en-dbp}); draw(); }
    else updateCross(); });
  const endDrag = () => { drag=null; viewEl.classList.remove("grabbing"); };
  viewEl.addEventListener("pointerup", endDrag); viewEl.addEventListener("pointercancel", endDrag);
  viewEl.addEventListener("pointerleave", () => { lastCX=null; crossEl.style.display="none"; });
  viewEl.addEventListener("dblclick", e => { e.preventDefault(); zoomAt(xToBp(e.clientX), e.shiftKey?1.8:0.55); });
  viewEl.addEventListener("keydown", e => { const sp=view.end-view.start, mid=(view.start+view.end)/2; let h=true;
    if(e.key==="ArrowLeft") animateTo({start:view.start-sp*0.15, end:view.end-sp*0.15});
    else if(e.key==="ArrowRight") animateTo({start:view.start+sp*0.15, end:view.end+sp*0.15});
    else if(e.key==="ArrowUp"||e.key==="+"||e.key==="=") zoomAt(mid, 0.625);
    else if(e.key==="ArrowDown"||e.key==="-"||e.key==="_") zoomAt(mid, 1.6);
    else if(e.key==="Home"||e.key==="0") animateTo({start:0, end:L});
    else if(e.key==="PageUp") animateTo({start:view.start-sp*0.5, end:view.end-sp*0.5});
    else if(e.key==="PageDown") animateTo({start:view.start+sp*0.5, end:view.end+sp*0.5});
    else h=false;
    if(h) e.preventDefault(); });
  host.querySelectorAll("[data-th]").forEach(b => b.onclick = () => { _figState[key+"_th"]=b.dataset.th; mountGenomeViewer(host, model, opts); });
  host.querySelector('[data-z="in"]').onclick  = () => zoomAt((view.start+view.end)/2, 0.625);
  host.querySelector('[data-z="out"]').onclick = () => zoomAt((view.start+view.end)/2, 1.6);
  host.querySelector('[data-z="fit"]').onclick = () => animateTo({start:0, end:L});
  // export what is on screen (WYSIWYG): the current zoom/pan window, not always the whole element (press Home/fit for the full view)
  host.querySelector('[data-exp="svg"]').onclick = () => downloadSVG(svgGenome(model, {start:view.start,end:view.end}, 920, "transparent", true), base);
  host.querySelector('[data-exp="png"]').onclick = () => downloadPNG(svgGenome(model, {start:view.start,end:view.end}, 920, "transparent", true), base);
  requestAnimationFrame(draw);
}

/* (svgStructure / svgGeneModel replaced by the interactive genome-browser viewer above) */

// gel render palettes — transparent (default), dark, white, UV transilluminator, monochrome
const GELPAL = {
  transparent: { paper:"none",    gel:"#0f1316", well:"#04060a", stroke:"#2a3138", ink:"#5a656f", on:OK.on,     off:OK.off,    ladder:OK.ladder, glow:1.4, band:2.6 },
  dark:        { paper:"#0b0e11", gel:"#0f1316", well:"#04060a", stroke:"#232a30", ink:"#8792a0", on:OK.on,     off:OK.off,    ladder:OK.ladder, glow:1.4, band:2.6 },
  white:       { paper:"#ffffff", gel:"#ededed", well:"#c4c4c4", stroke:"#cccccc", ink:"#555555", on:"#151515", off:"#992222", ladder:"#9a9a9a", glow:0.3, band:2.6 },
  uv:          { paper:"#050310", gel:"#0a0714", well:"#000000", stroke:"#1c1236", ink:"#9fb4d8", on:"#5bff6b", off:"#ffcf47", ladder:"#79d0ff", glow:3.2, band:3.1 },
  mono:        { paper:"#0d0d0d", gel:"#181818", well:"#000000", stroke:"#2b2b2b", ink:"#b2b2b2", on:"#f2f2f2", off:"#9a9a9a", ladder:"#cfcfcf", glow:2.0, band:2.9 },
};
function svgGel(data, bg){
  // data: { lanes:[{label, amplicons}] } or legacy { amplicons }
  const lanes = data.lanes || [{ label: data.laneLabel || "PCR", amplicons: (data.amplicons||[]) }];
  const sizes = lanes.flatMap(l => (l.amplicons||[]).map(a=>a.length));
  const smallest = sizes.length ? Math.min(...sizes) : 90;
  const minbp = Math.max(25, Math.min(90, smallest - 10)), maxbp = Math.max(1600, ...sizes);   // floor tracks the smallest band so 70-89 bp bands resolve
  const LADDER = [1500,1000,700,500,400,300,200,100,50].filter(m => m >= minbp && m <= maxbp);
  const laneW = 40, gap = 12, x0 = 62, top = 48, botPad = 46, H = 366, bot = H - botPad;
  const y = bp => top + (Math.log(maxbp)-Math.log(Math.max(bp,minbp)))/(Math.log(maxbp)-Math.log(minbp))*(bot-top);
  const cols = 1 + lanes.length;                              // ladder lane + sample lanes
  const W = Math.max(x0 + cols*(laneW+gap) + 12, 300);        // min width so the legend never clips
  const P = GELPAL[bg] || GELPAL.transparent;
  const laneX = i => x0 + i*(laneW+gap);
  let s = `<svg xmlns="http://www.w3.org/2000/svg" viewBox="0 0 ${W} ${H}" width="${W}" height="${H}" font-family="${FIGFONT}">`;
  s += `<defs><filter id="glow" x="-40%" y="-140%" width="180%" height="380%"><feGaussianBlur stdDeviation="${P.glow}" result="b"/><feMerge><feMergeNode in="b"/><feMergeNode in="SourceGraphic"/></feMerge></filter></defs>`;
  if(P.paper!=="none") s += `<rect width="${W}" height="${H}" fill="${P.paper}"/>`;
  // gel body (spans all lanes)
  s += `<rect x="${(x0-7).toFixed(1)}" y="${top-16}" width="${(cols*(laneW+gap)+2).toFixed(1)}" height="${(bot-top+30).toFixed(1)}" rx="2" fill="${P.gel}" stroke="${P.stroke}"/>`;
  s += `<text x="${x0-14}" y="${top-19}" fill="${P.ink}" font-size="8" text-anchor="end">bp</text>`;
  // MW labels aligned to ladder rungs
  LADDER.forEach(m => { const yy=y(m); s += `<text x="${x0-14}" y="${(yy+2.5).toFixed(1)}" fill="${P.ink}" font-size="8" text-anchor="end">${m}</text>`; });
  const drawLane = (col, label, bands, isLadder) => {
    const lx = laneX(col);
    s += `<rect x="${(lx+4).toFixed(1)}" y="${top-13}" width="${laneW-8}" height="4" rx="1" fill="${P.well}"/>`;   // well slot (to scale)
    s += `<text x="${(lx+laneW/2).toFixed(1)}" y="${top-18}" fill="${P.ink}" font-size="8.5" text-anchor="middle">${esc(label)}</text>`;
    (bands||[]).forEach(bd => { const yy=y(bd.size), col2=bd.color, h=isLadder?1.6:P.band;
      s += `<rect x="${(lx+3).toFixed(1)}" y="${(yy-h/2).toFixed(1)}" width="${laneW-6}" height="${h}" rx="1" fill="${col2}"${isLadder?"":' filter="url(#glow)"'}><title>${bd.size} bp${bd.t?" · "+bd.t:""}</title></rect>`; });
    if(!isLadder && !(bands||[]).length) s += `<text x="${(lx+laneW/2).toFixed(1)}" y="${bot+13}" fill="${P.ink}" font-size="7" text-anchor="middle">—</text>`;
  };
  drawLane(0, "L", LADDER.map(m=>({size:m, color:P.ladder})), true);
  lanes.forEach((l,i) => drawLane(i+1, l.label,
    (l.amplicons||[]).map(a=>({size:a.length, color:a.on_target?P.on:P.off, t:(a.on_target?"on-target":"off-target")+" · "+esc(a.source)})), false));
  // self-describing legend, well clear of the gel body
  const ly = H - 12;
  s += `<circle cx="${x0}" cy="${ly}" r="3" fill="${P.on}"/><text x="${x0+7}" y="${ly+3}" fill="${P.ink}" font-size="8">on-target</text>`
    +  `<circle cx="${x0+64}" cy="${ly}" r="3" fill="${P.off}"/><text x="${x0+71}" y="${ly+3}" fill="${P.ink}" font-size="8">off-target</text>`
    +  `<text x="${x0+140}" y="${ly+3}" fill="${P.ink}" font-size="8">L = MW ladder (bp)</text>`;
  return s + "</svg>";
}

/* interactive figure viewer: transparent by default, zoom (wheel/buttons) + pan (drag) + fit, resizable, first-class export */
const _figView = {};                                  // key -> {scale, tx, ty}
const _MODELABEL = { dark:"dark", white:"white", uv:"UV", mono:"mono", transparent:"transparent" };
function _svgWH(svg){ const m = svg.match(/viewBox="0 0 ([\d.]+) ([\d.]+)"/); return { w:+(m?m[1]:800), h:+(m?m[2]:400) }; }
function mountFigure(host, buildFn, baseName, opts){
  if(!host) return;
  delete _figView[host.id || baseName];               // fresh render → refit to the new figure
  _renderFig(host, buildFn, baseName, opts || {});
}
function _renderFig(host, buildFn, baseName, opts){
  const key = host.id || baseName;
  const modes = opts.modes || ["dark","white"];
  const dflt = modes[0];                                // themed view default (e.g. a dark gel), never a checkerboard
  const bg = _figState[key] || dflt;
  const svg = buildFn(bg);
  const viewH = opts.viewH || 360;
  const modeBtns = modes.map(m => `<button class="btn sm" data-bg="${m}">${_MODELABEL[m]||m}</button>`).join("");
  host.innerHTML = `<div class="figtoolbar">
      <span class="lbl">bg</span>${modeBtns}
      <span class="figmode" data-mode>view: ${_MODELABEL[bg]||bg} · export: transparent</span>
      <span class="grow"></span>
      <button class="btn sm" data-z="out" title="zoom out" aria-label="zoom out">−</button>
      <button class="btn sm" data-z="fit" title="fit to view" aria-label="fit figure to view">fit</button>
      <button class="btn sm" data-z="in" title="zoom in" aria-label="zoom in">+</button>
      <button class="btn sm" data-exp="svg" title="download vector SVG">⭳ SVG</button>
      <button class="btn sm" data-exp="png" title="download high-res PNG">⭳ PNG</button>
    </div>
    <div class="figview${opts.cls?' '+opts.cls:''}${bg==='transparent'?' clearbg':''}" data-view style="height:${viewH}px">
      <div class="figpan" data-pan>${svg}</div>
    </div>
    <div class="fighint">scroll = zoom · drag = pan · <b>fit</b> resets · drag the panel's bottom edge to resize</div>`;
  const view = host.querySelector("[data-view]"), pan = host.querySelector("[data-pan]");
  const size = _svgWH(svg);
  const apply = () => { const v=_figView[key]; if(v) pan.style.transform = `translate(${v.tx.toFixed(1)}px,${v.ty.toFixed(1)}px) scale(${v.scale.toFixed(4)})`; };
  const fit = () => { const vw=view.clientWidth||600, vh=view.clientHeight||viewH;
    const scale = Math.min(vw/size.w, vh/size.h)*0.94;
    _figView[key] = { scale, tx:(vw-size.w*scale)/2, ty:(vh-size.h*scale)/2 }; apply(); };
  const zoomAt = (cx,cy,f) => { const v=_figView[key]||{scale:1,tx:0,ty:0}; const ns=Math.min(16,Math.max(0.08,v.scale*f));
    v.tx = cx-(cx-v.tx)*(ns/v.scale); v.ty = cy-(cy-v.ty)*(ns/v.scale); v.scale=ns; _figView[key]=v; apply(); };
  host.querySelectorAll("[data-bg]").forEach(b => { if(b.dataset.bg===bg) b.classList.add("primary");
    b.onclick = () => { _figState[key] = b.dataset.bg; _renderFig(host, buildFn, baseName, opts); }; });
  host.querySelector('[data-z="in"]').onclick  = () => { const r=view.getBoundingClientRect(); zoomAt(r.width/2, r.height/2, 1.25); };
  host.querySelector('[data-z="out"]').onclick = () => { const r=view.getBoundingClientRect(); zoomAt(r.width/2, r.height/2, 0.8); };
  host.querySelector('[data-z="fit"]').onclick = fit;
  host.querySelector('[data-exp="svg"]').onclick = () => downloadSVG(buildFn("transparent"), baseName);   // export transparent for publication
  host.querySelector('[data-exp="png"]').onclick = () => downloadPNG(buildFn("transparent"), baseName);
  view.addEventListener("wheel", e => { e.preventDefault(); const r=view.getBoundingClientRect(); zoomAt(e.clientX-r.left, e.clientY-r.top, e.deltaY<0?1.12:0.89); }, {passive:false});
  let drag=null;
  view.addEventListener("pointerdown", e => { const v=_figView[key]||{tx:0,ty:0}; drag={x:e.clientX,y:e.clientY,tx:v.tx,ty:v.ty}; view.setPointerCapture(e.pointerId); view.classList.add("grabbing"); });
  view.addEventListener("pointermove", e => { if(!drag) return; const v=_figView[key]; v.tx=drag.tx+(e.clientX-drag.x); v.ty=drag.ty+(e.clientY-drag.y); apply(); });
  const end = () => { drag=null; view.classList.remove("grabbing"); };
  view.addEventListener("pointerup", end); view.addEventListener("pointercancel", end);
  if(_figView[key]) apply(); else requestAnimationFrame(fit);   // fresh render fits after layout; mode toggle preserves view
}
function downloadBlob(blob, name){ const u=URL.createObjectURL(blob), a=document.createElement("a"); a.href=u; a.download=name; a.click(); setTimeout(()=>URL.revokeObjectURL(u), 2000); }
/* ---------- tabular exports: CSV (Excel-ready) and TSV for every result table ---------- */
function csvEscape(v, sep){ v = String(v == null ? "" : v);
  if(/^[=+\-@\t\r]/.test(v)) v = "'" + v;                    // neutralize spreadsheet formula injection (CWE-1236)
  return (v.includes(sep) || v.includes('"') || v.includes("\n")) ? '"' + v.replace(/"/g, '""') + '"' : v; }
function downloadTable(headers, rows, base, fmt){
  const sep = fmt === "tsv" ? "\t" : ",";
  const lines = [headers.map(h => csvEscape(h, sep)).join(sep)].concat(rows.map(r => r.map(c => csvEscape(c, sep)).join(sep)));
  downloadBlob(new Blob(["﻿" + lines.join("\r\n")], { type: fmt === "tsv" ? "text/tab-separated-values" : "text/csv" }), base + "." + (fmt === "tsv" ? "tsv" : "csv"));
}
window.exportTable = (kind, fmt) => {
  let headers, rows, base = "TEagle";
  if(kind === "struct"){ const s = (state.lastRec && state.lastRec.structural) || [];
    headers = ["type","start","end","span_length","repeat_or_feature_length","identity_or_metric","method"];
    rows = s.map(e => { const sp = e.element_span || e.five_prime || e.pos || e.upstream || [null,null];
      const arm = e.ltr_len || e.tir_len || e.length || "";
      const span = (sp[0]!=null && sp[1]!=null) ? sp[1]-sp[0] : "";
      return [e.type, sp[0], sp[1], span, arm, e.identity!=null?e.identity+"%":(e.motif||""), e.method||""]; }); base = "TEagle_structural"; }
  else if(kind === "orf"){ const o = state.orfs || []; headers = ["strand","frame","start","end","length_aa"];
    rows = o.map(x => [x.strand, x.frame, x.start, x.end, x.length_aa]); base = "TEagle_orfs"; }
  else if(kind === "domain"){ const d = state.domains || []; headers = ["domain","label","pfam","aa_start","aa_end","nt_start","nt_end","score","evalue"];
    rows = d.map(x => [x.domain, x.label, x.pfam||"", x.aa[0], x.aa[1], x.nt[0], x.nt[1], x.score, x.evalue]); base = "TEagle_domains"; }
  else if(kind === "family"){ const f = state.family || []; headers = ["class_family","family","q_start","q_end","strand","divergence_pct","score"];
    rows = f.map(x => [x.class_family, x.family, x.q_start, x.q_end, x.strand, x.divergence, x.score]); base = "TEagle_dfam_families"; }
  else if(kind === "primers"){ const c = state.candidates || []; headers = ["id","fwd_5to3","rev_5to3","product_size_bp","left_tm","right_tm","left_gc_pct","right_gc_pct","penalty"];
    rows = c.map(x => [x.id, x.left_seq, x.right_seq, x.product_size, x.left_tm, x.right_tm, x.left_gc, x.right_gc, x.penalty]); base = "TEagle_primers"; }
  else if(kind === "pcr"){ const a = (state.lastPcr && state.lastPcr.amplicons) || []; headers = ["pair","source","start","end","length_bp","fwd_mismatches","rev_mismatches","call"];
    rows = a.map(x => [x.pair||"", x.source, x.start, x.end, x.length, x.fwd_mm, x.rev_mm, x.on_target?"on-target":"off-target"]); base = "TEagle_pcr_amplicons"; }
  else if(kind === "gene"){ const g = state.features; if(!g) return; const seg = [];
    (g.exons||[]).forEach(e => seg.push(["exon", e.start, e.end, e.end-e.start, e.strand]));
    (g.introns||[]).forEach(x => seg.push(["intron", x.start, x.end, x.end-x.start, x.strand]));
    (g.cds||[]).forEach(c => seg.push(["CDS", c.start, c.end, c.end-c.start, c.strand]));
    seg.sort((a,b) => a[1]-b[1]); headers = ["segment","start","end","length","strand"]; rows = seg; base = "TEagle_gene_structure"; }
  else if(kind === "splice"){ const sp = state.splice; if(!sp) return; headers = ["intron","start","end","length","donor","acceptor","canonical"];
    rows = (sp.gm.introns||[]).map((i,k) => [k+1, i.start, i.end, i.end-i.start, i.donor, i.acceptor, i.canonical]); base = "TEagle_splice_introns"; }
  else return;
  if(!rows.length) return;
  downloadTable(headers, rows, base, fmt);
};
function csvBtns(kind){ return `<button class="btn sm" onclick="exportTable('${kind}','csv')" title="download as CSV (opens in Excel)">⭳ CSV</button><button class="btn sm" onclick="exportTable('${kind}','tsv')" title="download as tab-separated values">⭳ TSV</button>`; }
function downloadSVG(svg, base){ downloadBlob(new Blob([svg], {type:"image/svg+xml"}), base+".svg"); }
function downloadPNG(svg, base, scale=3){
  const vb = svg.match(/viewBox="0 0 (\d+) (\d+)"/); const W=+(vb?vb[1]:800), H=+(vb?vb[2]:400);
  const img = new Image();
  img.onload = () => { const cv=document.createElement("canvas"); cv.width=W*scale; cv.height=H*scale;
    const x=cv.getContext("2d"); x.scale(scale,scale); x.drawImage(img,0,0); cv.toBlob(b=>downloadBlob(b, base+".png"), "image/png"); };
  img.src = "data:image/svg+xml;base64," + btoa(unescape(encodeURIComponent(svg)));
}

/* ---------- collapsible cards: all result cards collapsed by default, revealed as analysis progresses ---------- */
document.querySelectorAll(".main .card").forEach(card => {
  const ch = card.querySelector(".ch");
  if(ch && !ch.querySelector(".caret")){
    const car = document.createElement("span"); car.className = "caret"; car.textContent = "▾"; ch.appendChild(car);
    ch.setAttribute("role", "button"); ch.tabIndex = 0; ch.setAttribute("aria-expanded", "false");
    const toggle = () => { const c = card.classList.toggle("collapsed"); ch.setAttribute("aria-expanded", String(!c)); };
    ch.addEventListener("click", toggle);
    ch.addEventListener("keydown", e => { if(e.key==="Enter"||e.key===" "){ e.preventDefault(); toggle(); } });
  }
  card.classList.add("collapsed");
});
(function(){                                   // environment panel (rail): collapsible, collapsed by default
  const kv = document.getElementById("envkv"), hdr = kv && kv.previousElementSibling;
  if(hdr){
    hdr.classList.add("envhdr", "closed"); kv.classList.add("collapsed");
    const car = document.createElement("span"); car.className = "caret"; car.textContent = "▾"; hdr.appendChild(car);
    hdr.setAttribute("role", "button"); hdr.tabIndex = 0; hdr.setAttribute("aria-expanded", "false");
    const t = () => { const c = kv.classList.toggle("collapsed"); hdr.classList.toggle("closed"); hdr.setAttribute("aria-expanded", String(!c)); };
    hdr.addEventListener("click", t);
    hdr.addEventListener("keydown", e => { if(e.key==="Enter"||e.key===" "){ e.preventDefault(); t(); } });
  }
})();
/* accessibility: label controls + live region for errors/warnings */
document.querySelectorAll("input, textarea").forEach(el => { if(!el.getAttribute("aria-label") && el.placeholder) el.setAttribute("aria-label", el.placeholder); });
(function(){ const eb = document.getElementById("errbox"); if(eb){ eb.setAttribute("role","alert"); eb.setAttribute("aria-live","assertive"); } })();
document.querySelectorAll(".pf").forEach(pf => {
  const l = pf.querySelector("label"), i = pf.querySelector("input, select");
  if(l && i){ if(!i.id) i.id = "pf_" + l.textContent.trim().replace(/\W+/g,"_"); l.setAttribute("for", i.id); i.setAttribute("aria-label", l.textContent.trim()); }
});
function openMenu(x,y,title,items,trigger){
  closeMenu();
  _menuTrigger = trigger || null;
  _menu = document.createElement("div"); _menu.className = "ctxmenu"; _menu.setAttribute("role","menu");
  _menu.innerHTML = `<div class="hd">${title}</div>` + items.map((it,i)=>`<button role="menuitem" data-i="${i}">${it.label}</button>`).join("");
  document.body.appendChild(_menu);
  _menu.style.left = Math.min(Math.max(x,4), innerWidth - _menu.offsetWidth - 8) + "px";
  _menu.style.top = Math.min(Math.max(y,4), innerHeight - _menu.offsetHeight - 8) + "px";
  const btns = [..._menu.querySelectorAll("button")];
  btns.forEach((b,i)=> b.onclick = ev => { ev.stopPropagation(); items[i].fn(b); setTimeout(()=>closeMenu(), 1000); });  // dismiss after activation (toast shows first)
  _menu.addEventListener("keydown", e => {
    const idx = btns.indexOf(document.activeElement);
    if(e.key==="ArrowDown"){ e.preventDefault(); btns[(idx+1)%btns.length].focus(); }
    else if(e.key==="ArrowUp"){ e.preventDefault(); btns[(idx-1+btns.length)%btns.length].focus(); }
    else if(e.key==="Escape" || e.key==="Tab"){ e.preventDefault(); e.stopPropagation(); closeMenu(true); }  // Tab must not orphan the menu
  });
  if(btns[0]) btns[0].focus();                          // move focus into the menu (keyboard-operable)
}
window.featMenu = (ev, start, end, strand, label, extra) => {
  ev.stopPropagation();
  const dna = strand === "-" ? revcomp(slice(start,end)) : slice(start,end);
  const items = [
    {label:"⧉  Copy FASTA", fn:b=>copyText(`>${label}_${start}-${end}${strand==="-"?"_rev":""}\n${dna}`, b)},
    {label:"⧉  Copy DNA", fn:b=>copyText(dna, b)},
    {label:`⧉  Copy coords (${start}–${end} ${strand})`, fn:b=>copyText(`${start}-${end} ${strand}`, b)},
    {label:"⌖  Design primer here", fn:()=>{ closeMenu(); designForDomain(start,end,label); }},
  ];
  if(extra && extra.protein) items.splice(2, 0, {label:"⧉  Copy protein", fn:b=>copyText(extra.protein, b)});
  let x = ev.clientX, y = ev.clientY;
  const trig = ev.currentTarget && ev.currentTarget.getBoundingClientRect ? ev.currentTarget : document.activeElement;
  if((!x && !y) && trig && trig.getBoundingClientRect){ const r = trig.getBoundingClientRect(); x = r.left + 24; y = r.bottom - 4; }  // keyboard: anchor to row
  openMenu(x, y, label, items, trig);
};
window.copyAllFasta = (kind) => {
  let recs = [];
  if(kind === "orf") recs = (state.orfs||[]).map(o => `>ORF_${o.strand}${o.frame}_${o.start}-${o.end}_${o.length_aa}aa\n${o.strand==="-"?revcomp(slice(o.start,o.end)):slice(o.start,o.end)}`);
  else if(kind === "struct") recs = (state.structural||[]).map(e => { const sp=e.element_span||e.five_prime||e.pos||[0,0]; return `>${e.type.split(" ")[0]}_${sp[0]}-${sp[1]}\n${slice(sp[0],sp[1])}`; });
  else if(kind === "family") recs = (state.family||[]).map(h => `>${h.family}_${h.class_family.replace("/","-")}_${h.q_start}-${h.q_end}\n${h.strand==="-"?revcomp(slice(h.q_start,h.q_end)):slice(h.q_start,h.q_end)}`);
  copyText(recs.join("\n"), event && event.currentTarget);
};
window.copyAmplicons = () => {
  const a = (state.lastPcr && state.lastPcr.amplicons) || [];
  copyText(a.map(x => `>amplicon_${String(x.source).replace(/\W+/g,"_")}_${x.start}-${x.end}_${x.length}bp_${x.on_target?"on":"off"}target\n${x.seq||""}`).join("\n"), event && event.currentTarget);
};

/* ---------- health ---------- */
(async () => {
  try {
    const h = await (await fetch("/api/health")).json();
    $("#status").classList.add("live");
    $("#statusTxt").textContent = `backend live · primer3 ${h.primer3}`;
    if (h.core) $("#ver").textContent = "v" + h.core;   // single source of truth: teagle_core.__version__
  } catch { $("#statusTxt").textContent = "backend offline"; }
})();

/* ---------- environment status ---------- */
(async () => {
  try {
    const e = await (await fetch("/api/env")).json();
    const pyOk = e.python_ok ? "ok" : "old";
    const pkgs = (e.packages||[]).map(p => `${p.name} ${p.ok?"✓":"✗ "+(p.installed||"missing")}`);
    const state = e.error ? `<span style="color:var(--bad)">${e.error.slice(0,40)}</span>`
      : e.needs_install ? `<span style="color:var(--warn)">install needed</span>`
      : `<span style="color:var(--good)">up to date</span>`;
    const bw = e.backends||{};
    $("#envkv").innerHTML = `
      <dt>state</dt><dd>${state}${e.first_run?" · first run":""}</dd>
      <dt>python</dt><dd>${e.python} <span style="color:${e.python_ok?'var(--good)':'var(--bad)'}">${pyOk}</span></dd>
      <dt>packages</dt><dd style="font-size:10.5px">${pkgs.join("<br>")}</dd>
      <dt>webview2</dt><dd>${bw.webview2?'<span style="color:var(--good)">present</span>':'<span style="color:var(--faint)">—</span>'}</dd>
      <dt>wsl2</dt><dd>${bw.wsl2==="available"?'<span style="color:var(--good)">available</span>':`<span style="color:var(--faint)">${bw.wsl2||"—"}</span>`}</dd>
      <dt>signature</dt><dd class="hash">${e.signature||"—"}</dd>`;
  } catch { $("#envkv").innerHTML = '<dt>status</dt><dd style="color:var(--bad)">env check unavailable</dd>'; }
})();

/* WSL2 is optional — only Dfam naming / de-novo splice need it; core classification does not. */
function wslAbsentHtml(feature){
  return `<div class="wsl-absent"><b>WSL2 is not installed</b> — this step is <b>optional</b>. `
    + `WSL2 is Windows' built-in Linux layer, needed only for ${feature}; the domain-based superfamily `
    + `classification above works without it.<br>To enable it, open <b>PowerShell as Administrator</b>, run `
    + `<code>wsl --install</code> <button class="btn sm" onclick="navigator.clipboard&&navigator.clipboard.writeText('wsl --install')" title="copy command">⧉ copy</button>, `
    + `then restart Windows. <a href="https://learn.microsoft.com/windows/wsl/install" target="_blank" rel="noopener">Microsoft's WSL guide ↗</a></div>`;
}

/* ---------- WSL family-annotation backend (Dfam / RepeatMasker) ---------- */
async function initWsl(){
  const st = $("#wslStatus");
  try {
    const w = await (await fetch("/api/wsl/status")).json();
    if(w.error){ st.innerHTML = `<span style="color:var(--bad)">WSL status error: ${w.error}</span>`; return; }
    if(!w.wsl2){ st.innerHTML = wslAbsentHtml("Dfam family-level naming"); return; }
    if(w.ready){
      st.innerHTML = `<span style="color:var(--good)">● ready</span> · RepeatMasker ${w.repeatmasker} · Dfam curated installed · distro ${w.distro}`;
      $("#annotate").disabled = false;
    } else {
      st.innerHTML = `WSL2 ok (${w.distro}); annotation stack not installed (RepeatMasker ${w.repeatmasker||"missing"}, Dfam ${w.dfam?"ok":"missing"}). ` +
        `<button class="btn sm" id="wslInstall" style="margin-left:6px">Install backend</button>`;
      const ib = $("#wslInstall"); if(ib) ib.onclick = startInstall;
    }
  } catch(e){ st.innerHTML = `<span style="color:var(--bad)">WSL check failed: ${e.message}</span>`; }
}
async function startInstall(){
  const eta = await etaText("install", 0);
  $("#wslStatus").innerHTML = `installing WSL backend (RepeatMasker + Dfam curated + minimap2) — ${eta||"this can take several minutes"}…`;
  try { await api("/api/wsl/install", {}); } catch(e){}
  const poll = setInterval(async () => {
    try {
      const r = await (await fetch("/api/wsl/install_log")).json();
      const log = r.log || "";
      $("#wslStatus").innerHTML = `<pre style="font-family:var(--mono);font-size:10px;white-space:pre-wrap;color:var(--dim);max-height:120px;overflow:auto">${log.slice(-700)}</pre>`;
      if(log.includes("[teagle] DONE")){ clearInterval(poll); initWsl(); }
      else if(log.includes("[teagle] FAILED")){ clearInterval(poll);
        const m = (log.match(/\[teagle\] FAILED:[^\n]*/) || ["install failed"])[0];
        $("#wslStatus").innerHTML = `<span style="color:var(--bad)">✕ ${m}</span> <button class="btn sm" id="wslRetry">Retry install</button>`;
        const rb = $("#wslRetry"); if(rb) rb.onclick = startInstall;
      }
    } catch(e){}
  }, 5000);
}
/* panel-03 sequence source: loaded specimen | upload FASTA | paste (independent of panel 01) */
function wslResolveSeq(){
  const src = $("#wslSource") ? $("#wslSource").value : "specimen";
  if(src === "paste")  return { seq: ($("#wslPaste").value||""), label: "pasted sequence", useSource: false };
  if(src === "upload") return { seq: (state.wslSeq||""), label: state.wslSeqName || "uploaded file", useSource: false };
  return { seq: (state.seq || seqEl.value || ""), label: "loaded specimen", useSource: true };
}
if($("#wslSource")) $("#wslSource").onchange = e => {
  const v = e.target.value, meta = $("#wslSrcMeta");
  $("#wslPaste").hidden = v !== "paste";
  if(v === "upload"){ $("#wslFile").click(); }
  else if(v === "specimen"){ meta.textContent = state.seq ? "using the sequence loaded in panel 01" : "no specimen loaded yet in panel 01"; }
  else if(v === "paste"){ meta.textContent = "paste a FASTA or raw sequence below"; $("#wslPaste").focus(); }
};
if($("#wslFile")) $("#wslFile").onchange = async (e) => {
  const f = e.target.files && e.target.files[0], meta = $("#wslSrcMeta");
  if(!f){ if(!state.wslSeq) $("#wslSource").value = "specimen"; return; }
  meta.textContent = `reading ${f.name}…`;
  try {
    let text;
    if(f.name.toLowerCase().endsWith(".gz")){
      if(typeof DecompressionStream === "undefined") throw new Error("gzip not supported in this browser");
      text = await new Response(f.stream().pipeThrough(new DecompressionStream("gzip"))).text();
    } else { text = await f.text(); }
    if(!text.trim()) throw new Error("file is empty");
    state.wslSeq = text; state.wslSeqName = f.name;
    const nrec = (text.match(/^>/gm)||[]).length, bp = text.replace(/^>.*$/gm,"").replace(/\s/g,"").length;
    meta.innerHTML = `<b>${f.name}</b> · ${bp.toLocaleString()} bp${nrec>1?` · ${nrec} records (first record annotated)`:""} — ready to annotate`;
  } catch(er){ meta.textContent = "✕ upload failed — " + er.message; state.wslSeq = ""; $("#wslSource").value = "specimen"; }
  finally { e.target.value = ""; }
};
if($("#annotate")) $("#annotate").onclick = async () => {
  clearErr();
  const sp = $("#wslSpecies").value.trim();
  const r = wslResolveSeq();
  if(!r.seq.trim()){ err(`No sequence for family annotation — the ${r.label} is empty. Choose a source above.`); return; }
  const nrec = (r.seq.match(/^>/gm)||[]).length;
  const bp = r.seq.replace(/^>.*$/gm,"").replace(/\s/g,"").length;
  const btn = $("#annotate"), o = btn.textContent; btn.disabled = true; btn.classList.add("pending"); btn.textContent = "◴ running…";
  const eta = await etaText("annotate", bp);
  $("#wslBody").innerHTML = `<div class="mini">running RepeatMasker against Dfam on the ${r.label}${nrec>1?` (first of ${nrec} records)`:""}${sp?` (species: ${sp})`:""}…<br>${eta?eta+" · ":""}<span class="etaclock">0s elapsed</span></div>`;
  const stop = liveTimer($("#wslBody"));
  try {
    const d = await api("/api/annotate", { sequence: r.seq, species: sp, source: r.useSource ? state.source : null, timeout: 600 });
    if(!d.ok){ $("#wslBody").innerHTML = `<div class="err-banner">${d.error}</div>`; return; }
    renderFamily(d);
    if(d.elapsed_s) $("#wslMeta").textContent += ` · ${d.elapsed_s}s`;
    if(d.provenance) renderProvenance(d.provenance);
  } catch(e){ err("Family annotation failed — " + e.message); }
  finally { stop(); btn.disabled = false; btn.classList.remove("pending"); btn.textContent = o; }
};
const _LC = new Set(["Low_complexity", "Simple_repeat", "Satellite", "Unknown", "Unspecified"]);
function famTag(cf){
  const c = (cf.split("/")[0] || "");
  const cls = c === "LTR" ? "t-ltr" : c === "DNA" ? "t-tir" : c === "LINE" ? "t-orf"
            : c === "SINE" ? "t-good" : "t-plain";
  return `<span class="tag ${cls}"><span class="d" style="background:currentColor"></span>${cf}</span>`;
}
function renderFamily(d){
  expandCard($("#wslBody").closest(".card"));
  const te = d.hits.filter(h => !_LC.has(h.class_family));   // TE families only, not low-complexity
  const lc = d.hits.length - te.length;
  const noSpecies = !d.species || d.species.startsWith("(all");
  if(!te.length){
    const hint = noSpecies
      ? " Set the organism/species above — RepeatMasker needs a lineage, and without one only low-complexity is reported."
      : " The family for this organism may need an additional Dfam taxon partition (curated covers well-studied lineages).";
    $("#wslBody").innerHTML = `<div class="empty">No TE family named under the current criteria`+
      `${lc?` (${lc} low-complexity / simple-repeat region${lc>1?"s":""} found)`:""}.${hint}</div>`;
    $("#wslMeta").textContent = `0 TE families${lc?` · ${lc} low-cplx`:""}`;
    return;
  }
  state.family = te;
  const fams = [...new Set(te.map(h=>h.family))];
  const rows = te.map((h,i)=>`<tr class="clickable" tabindex="0" role="button" onclick="featMenu(event, ${h.q_start}, ${h.q_end}, '${h.strand}', '${h.family}')">
    <td class="mono">${i+1}</td><td>${famTag(h.class_family)}</td><td class="mono">${h.family}</td>
    <td class="coord">${h.q_start.toLocaleString()}–${h.q_end.toLocaleString()}</td>
    <td class="mono">${h.strand}</td><td class="coord">${h.divergence}%</td><td class="coord">${h.score}</td></tr>`).join("");
  $("#wslBody").innerHTML = `
    <div class="classbn"><span class="big">Dfam · ${te[0].class_family}</span>
      <span class="cf cf-High">Detected</span>
      <span class="kls">${fams.join(" · ")} · Dfam 4.0 curated ${srcChip("Dfam")} · RepeatMasker ${d.repeatmasker_version} ${srcChip("RepeatMasker")}</span>
      <span class="exp">Family-level annotation from RepeatMasker (RMBLAST) against the Dfam curated library in WSL — species: ${d.species}. Layer-A homology (Detected), complementing the domain-based superfamily above.${lc?` (${lc} low-complexity region${lc>1?"s":""} omitted.)`:""}</span></div>
    <table><thead><tr><th>#</th><th>Class/family</th><th>Dfam family</th><th>Coords (0-based)</th><th>Str</th><th>Div</th><th>Score</th></tr></thead><tbody>${rows}</tbody></table>
    <div class="flexbtns" style="margin-top:8px"><button class="btn sm" onclick="copyAllFasta('family')">⧉ all family regions → FASTA</button>${csvBtns('family')}<span class="lbl">click any family row for its sequence</span></div>`;
  glossify($("#wslBody"));
  $("#wslMeta").textContent = `${te.length} family hits · RepeatMasker ${d.repeatmasker_version}`;
}
initWsl();

/* ---------- splice detection (de novo, minimap2 / WSL) ---------- */
async function initSplice(){
  const st = $("#spliceStatus"); if(!st) return;
  try {
    const w = await (await fetch("/api/wsl/status")).json();
    if(!w.wsl2){ st.innerHTML = wslAbsentHtml("de-novo splice detection"); return; }
    if(w.minimap2){ st.innerHTML = `<span style="color:var(--good)">● ready</span> · minimap2 ${w.minimap2} · align a transcript to resolve exon–intron structure`; $("#runSplice").disabled = false; $("#spliceHint").textContent = "paste a transcript, then detect"; }
    else { st.innerHTML = `minimap2 not installed in the WSL backend — it ships with the managed install (panel 03 “Install backend”).`; }
  } catch(e){ st.innerHTML = `<span style="color:var(--bad)">splice backend check failed: ${e.message}</span>`; }
}
if($("#runSplice")) $("#runSplice").onclick = async () => {
  clearErr();
  const genomic = state.seq || seqEl.value, tx = $("#spliceTx").value.trim();
  if(!genomic.trim()){ err("Load a genomic sequence first (fetch, upload, or paste, then Run analysis)."); return; }
  if(!tx){ err("Paste a transcript / cDNA / mRNA to align."); return; }
  const gbp = genomic.replace(/^>.*$/gm,"").replace(/\s/g,"").length;
  const btn = $("#runSplice"), o = btn.textContent; btn.disabled = true; btn.classList.add("pending"); btn.textContent = "◴ aligning…";
  const eta = await etaText("splice", gbp);
  $("#spliceBody").innerHTML = `<div class="mini">minimap2 -x splice aligning the transcript to the loaded sequence…<br>${eta?eta+" · ":""}<span class="etaclock">0s elapsed</span></div>`;
  const stop = liveTimer($("#spliceBody"));
  try {
    const d = await api("/api/splice", { sequence: genomic, transcript: tx, source: state.source, timeout: 300 });
    if(!d.ok){ $("#spliceBody").innerHTML = `<div class="err-banner">${d.error}</div>`; return; }
    renderSplice(d);
    if(d.provenance) renderProvenance(d.provenance);
  } catch(e){ err("Splice alignment failed — " + e.message); }
  finally { stop(); btn.disabled = false; btn.classList.remove("pending"); btn.textContent = o; }
};
function renderSplice(d){
  expandCard($("#spliceBody").closest(".card"));
  const len = (state.lastRec && state.lastRec.composition.length) || (seqEl.value||"").replace(/^>.*$/gm,"").replace(/\s/g,"").length || 1;
  const gm = { exons: d.exons, introns: d.introns, cds: [] };
  const rows = d.introns.map((i,k)=>`<tr class="clickable" tabindex="0" role="button" onclick="featMenu(event, ${i.start}, ${i.end}, '${d.strand}', 'intron_${k+1}')">
    <td class="mono">${k+1}</td><td class="coord">${i.start.toLocaleString()}–${i.end.toLocaleString()}</td><td class="coord">${(i.end-i.start).toLocaleString()}</td>
    <td class="mono">${i.donor}…${i.acceptor}</td>
    <td>${i.canonical?'<span class="tag t-good"><span class="d" style="background:var(--good)"></span>canonical</span>':'<span class="tag t-bad"><span class="d" style="background:var(--bad)"></span>non-canonical</span>'}</td></tr>`).join("")
    || `<tr><td colspan="5" class="mono" style="color:var(--faint);padding:10px">single exon — no introns detected</td></tr>`;
  $("#spliceBody").innerHTML = `
    <div class="classbn"><span class="big">${d.counts.exons} exon(s) · ${d.counts.introns} intron(s)</span>
      <span class="cf cf-High">de novo</span>
      <span class="kls">minimap2 -x splice · ${d.canonical_introns}/${d.counts.introns} canonical splice site(s) · strand ${d.strand} ${srcChip("minimap2")}</span>
      <span class="exp">Evidence-based exon–intron structure: the transcript was spliced-aligned to the loaded sequence; introns are alignment gaps, splice sites compared to canonical GT–AG / GC–AG / AT–AC motifs.</span></div>
    <div id="figSplice"></div>
    <table><thead><tr><th>#</th><th>Intron span (0-based)</th><th>Len</th><th>Splice site</th><th>Call</th></tr></thead><tbody>${rows}</tbody></table>
    <div class="flexbtns" style="margin-top:8px"><button class="btn sm" onclick="copySpliceFasta()">⧉ exons → FASTA</button>${csvBtns('splice')}
      <span class="lbl">exons/introns resolved by spliced alignment of your transcript</span></div>`;
  state.splice = { gm, len, strand: d.strand };
  glossify($("#spliceBody"));
  delete _gvState["figSplice"]; mountGenomeViewer($("#figSplice"), gvTracksFromGene(gm, len), { base:"TEagle_splice" });
}
window.copySpliceFasta = () => {
  const sp = state.splice; if(!sp) return;
  const minus = sp.strand === "-";
  copyText((sp.gm.exons||[]).map((x,i)=>{ const seq = minus ? revcomp(slice(x.start,x.end)) : slice(x.start,x.end);
    return `>exon_${i+1}_${x.start}-${x.end}_${sp.strand||"+"}\n${seq}`; }).join("\n"), event && event.currentTarget);
};
initSplice();

/* ---------- theme (persisted) ---------- */
(function(){
  const saved = localStorage.getItem("teagle_theme");
  if(saved) document.documentElement.setAttribute("data-theme", saved);
  else if(window.matchMedia && matchMedia("(prefers-color-scheme: light)").matches)
    document.documentElement.setAttribute("data-theme","light");
})();
$("#theme").onclick = () => {
  const r = document.documentElement;
  const next = r.getAttribute("data-theme")==="dark" ? "light" : "dark";
  r.setAttribute("data-theme", next);
  localStorage.setItem("teagle_theme", next);
};

/* ---------- specimen ---------- */
const seqEl = $("#seq");
function countNt(){
  const raw = seqEl.value.replace(/^>.*$/gm,"").replace(/\s/g,"");
  $("#charCount").textContent = raw.length.toLocaleString() + " nt";
}
seqEl.addEventListener("input", () => { countNt(); state.source = null; state.features = null; $("#accMeta").className = "accmeta"; $("#accMeta").innerHTML = "";
  if(state.analyzedSeq != null && seqEl.value !== state.analyzedSeq){          // edited after analysis -> downstream results are stale
    $("#design").disabled = true; $("#designHint").textContent = "sequence changed — re-run analysis";
    if($("#runpcr")) $("#runpcr").disabled = true;
  } });
function _staleBlock(){ if(state.analyzedSeq != null && seqEl.value !== state.analyzedSeq){
    err("The sequence changed since the last analysis. Click “Run analysis” again before designing primers or running PCR."); return true; } return false; }
$("#loadSample").onclick = () => { seqEl.value = mkSample(); countNt(); state.source = null; state.features = null; $("#accMeta").className="accmeta"; $("#accMeta").innerHTML=""; };
seqEl.value = mkSample(); countNt();

/* ---------- live accession fetch (real NCBI) ---------- */
async function doFetch(){
  const acc = $("#acc").value.trim();
  const meta = $("#accMeta");
  if(!acc){ meta.className="accmeta bad"; meta.textContent="enter an accession"; return; }
  meta.className = "accmeta loading"; meta.textContent = `fetching ${acc} from NCBI…`;
  $("#fetch").disabled = true;
  try{
    const d = await api("/api/fetch", { accession: acc });
    if(!d.ok){ meta.className="accmeta bad"; meta.textContent = "✕ " + d.error; return; }
    seqEl.value = d.fasta; countNt();
    state.source = { accession:d.accession, organism:d.organism, taxid:d.taxid,
                     source:d.source, endpoint:d.endpoint, retrievedUtc:d.retrievedUtc, fromCache:!!d.fromCache };
    state.features = d.features || null;
    if(d.organism){ const sp=$("#wslSpecies"); if(sp) sp.value = d.organism; }
    meta.className = "accmeta ok";
    const src = d.fromCache ? `<span title="reused from local cache — not re-downloaded">cached (local)</span>` : d.source||"NCBI";
    meta.innerHTML = `<b>${d.accession}</b> · ${d.organism||"?"} · taxid ${d.taxid||"?"} · ${src}<br>`+
      `${(d.title||"").slice(0,64)}<br><b>${(d.seq_length||d.length||0).toLocaleString()} bp</b> · ${d.moltype||"DNA"} · verify, then Run analysis`;
  }catch(e){ meta.className="accmeta bad"; meta.textContent = "✕ fetch failed — " + e.message; }
  finally{ $("#fetch").disabled = false; }
}
$("#fetch").onclick = doFetch;
$("#acc").addEventListener("keydown", e => { if(e.key==="Enter") doFetch(); });

/* ---------- upload FASTA (plain + gzip) ---------- */
if($("#uploadBtn")) $("#uploadBtn").onclick = () => $("#fileInput").click();
if($("#fileInput")) $("#fileInput").onchange = async (e) => {
  const f = e.target.files && e.target.files[0]; if(!f) return;
  const meta = $("#accMeta"); meta.className = "accmeta loading"; meta.textContent = `reading ${f.name}…`;
  try {
    let text;
    if(f.name.toLowerCase().endsWith(".gz")){
      if(typeof DecompressionStream === "undefined") throw new Error("gzip not supported in this browser");
      text = await new Response(f.stream().pipeThrough(new DecompressionStream("gzip"))).text();
    } else { text = await f.text(); }
    if(!text.trim()) throw new Error("file is empty");
    seqEl.value = text; countNt(); state.source = null; state.features = null;
    const nrec = (text.match(/^>/gm) || []).length;
    meta.className = "accmeta ok"; meta.innerHTML = `<b>${f.name}</b> loaded · ${(text.replace(/^>.*$/gm,"").replace(/\s/g,"").length).toLocaleString()} bp`+
      `${nrec>1?` · ${nrec} records (first is used downstream)`:""} · verify, then Run analysis`;
  } catch(er){ meta.className = "accmeta bad"; meta.textContent = "✕ upload failed — " + er.message; }
  finally { e.target.value = ""; }   // allow re-picking the same file
};

/* ---------- run analysis ---------- */
$("#run").onclick = async () => {
  clearErr();
  const btn = $("#run"); const orig = btn.textContent;
  btn.disabled = true; btn.classList.add("pending"); btn.textContent = "◴ analyzing…";
  try {
    state.seq = seqEl.value;
    const d = await api("/api/analyze", { sequence: state.seq, source: state.source });
    if(!d.records || !d.records.length){ err(d.warning || "Enter, upload, or fetch a sequence first."); return; }
    if(d.warning) $("#errbox").innerHTML = `<div class="pcrwarn">⚠ ${d.warning}</div>`;
    const rec = d.records[0];
    if(rec.notes && rec.notes.length) $("#errbox").innerHTML += `<div class="pcrwarn">⚠ ${rec.notes.join("; ")}</div>`;
    const c = rec.composition;
    $("#mLen").innerHTML = c.length.toLocaleString() + " <small>bp</small>";
    $("#mGC").innerHTML = c.gc + " <small>%</small>";
    $("#mN").innerHTML = c.n + " <small>%</small>";
    const vc = $("#mValidCell");
    vc.classList.toggle("ok", rec.valid); vc.classList.toggle("err", !rec.valid);
    $("#mValid").textContent = rec.valid ? "valid" : rec.invalid.length+" bad";
    const ncbi = state.source && state.source.accession
      ? srcChip("NCBI", "https://www.ncbi.nlm.nih.gov/nuccore/" + encodeURIComponent(state.source.accession))
      : (state.source ? "" : `<span class="srcnote" title="Sequence provided by you (upload or paste) — no external source">user-provided</span>`);
    $("#rRecords").innerHTML = `${d.records.length} ${ncbi}`;
    $("#rStruct").innerHTML = `${rec.structural.length} ${rec.structural.length ? srcChip("Wicker2007") : ""}`;
    $("#rOrf").textContent = rec.orfs.length;
    state.structural = rec.structural;
    const ltr = rec.structural.find(e => e.element_span);
    state.elemSpan = ltr ? ltr.element_span : null;
    renderStructure(rec);
    renderProvenance(d.provenance);
    $("#design").disabled = false; $("#designHint").textContent = "ready";
    state.analyzedSeq = state.seq;                     // snapshot for staleness detection
  } catch(e){ err("Analysis failed — " + e.message); }
  finally { btn.disabled = false; btn.classList.remove("pending"); btn.textContent = orig; }
};

/* ---------- structure render ---------- */
function tagFor(t){
  if(t.startsWith("LTR")) return `<span class="tag t-ltr"><span class="d" style="background:var(--ltr)"></span>LTR</span>`;
  if(t.startsWith("TIR")) return `<span class="tag t-tir"><span class="d" style="background:var(--tir)"></span>TIR</span>`;
  if(t.startsWith("TSD")) return `<span class="tag t-tsd"><span class="d" style="background:var(--tsd)"></span>TSD</span>`;
  if(t.startsWith("poly")) return `<span class="tag t-polya"><span class="d" style="background:var(--polya)"></span>tail</span>`;
  return `<span class="tag t-plain">feat</span>`;
}
const ORF_PP = 12;
function renderStructure(rec){
  const s = rec.structural;
  state.orfs = rec.orfs; state.orfPage = 0; state.domains = rec.domains || [];
  const cl = rec.classification, doms = rec.domains || [];
  $("#structMeta").textContent = (cl && cl.te_class !== "none" ? cl.te_class + " · " : "") +
    `${s.length} struct · ${doms.length} dom · ${rec.orfs.length} ORF`;
  const gm = (state.features && state.source) ? state.features : null;
  if(!(s.length || rec.orfs.length || doms.length || gm)){
    $("#structBody").innerHTML = `<div class="empty">No terminal repeats, domains, ORFs or tails detected in this sequence.</div>`;
    return;
  }
  const classBanner = (cl && cl.te_class !== "none") ? `<div class="classbn">
      <span class="big">${cl.te_class}</span>
      <span class="cf cf-${cl.confidence}">${cl.confidence} confidence</span>
      <span class="kls">${cl.class}${cl.order ? " · " + cl.order : ""}</span>
      <span class="exp">${cl.explanation}</span></div>` : "";
  const domPanel = doms.length ? `<div class="domarch">
      <div class="lbl" style="margin-bottom:8px;display:flex;align-items:center;gap:8px;flex-wrap:wrap">Domain architecture · ${doms.length} domains · 5′→3′ · <span style="color:var(--faint)">HMMER / Pfam — click a domain for its sequence</span>${csvBtns('domain')}</div>
      ${doms.map((d,i) => `<div class="domrow" tabindex="0" role="button" onclick="toggleDomain(${i})">
        <span class="dchip d-${d.domain}">${d.domain}</span>
        <span class="dl">${d.label}</span>
        <span class="dm">${d.pfam||""} · aa ${d.aa[0]}–${d.aa[1]} · nt ${d.nt[0]}–${d.nt[1]} · <span title="HMMER bit score — how strongly this region matches the domain profile; higher is stronger">score ${d.score}</span> · <span title="E-value: the number of matches this good expected by chance; lower is more significant (e.g. 1e-30 is highly significant)">E ${Number(d.evalue).toExponential(0)}</span> ${srcChip("Pfam")}</span>
        <span class="db"><button class="btn sm" onclick="event.stopPropagation();designForDomain(${d.nt[0]},${d.nt[1]},'${d.domain}')">⌖ domain PCR</button><span class="caret" id="dcar${i}">▸</span></span>
      </div><div class="domdetail" id="domdet${i}" style="display:none"></div>`).join("")}</div>` :
      `<div class="mini" style="margin-bottom:12px">No protein domains detected (no RT / integrase / transposase HMM hit above threshold — structural evidence only).</div>`;
  const rows = s.map(e => {
    const span = e.element_span || e.five_prime || e.pos || e.upstream || [null,null];
    const arm = e.ltr_len || e.tir_len;                                  // terminal-repeat arm length (LTR/TIR)
    let metric = e.identity!=null ? e.identity+"% id" : (e.length!=null ? e.length+" bp" : (e.motif? "motif "+e.motif : "—"));
    if(arm) metric = `${arm} bp repeat · ${metric}`;                     // disclose the arm here, so Len can mean the span consistently
    const len = (span[0]!=null && span[1]!=null) ? (span[1]-span[0]) : (e.length||"");   // Len == length of the displayed Span, for every row
    const click = span[0]!=null ? ` class="clickable" tabindex="0" role="button" onclick="featMenu(event, ${span[0]}, ${span[1]}, '+', '${e.type.split(" ")[0]}')"` : "";
    return `<tr${click}><td>${tagFor(e.type)}</td><td>${e.type}</td>
      <td class="coord">${span[0]}–${span[1]}</td><td class="coord">${len}</td>
      <td class="coord">${metric}</td><td class="mono" style="color:var(--faint);font-size:10.5px">${e.method||""} ${srcChip("Wicker2007")}</td></tr>`;
  }).join("") || `<tr><td colspan="6" class="mono" style="color:var(--faint);padding:10px">no terminal repeats or tails detected</td></tr>`;
  $("#structBody").innerHTML = classBanner + domPanel + `
    <div id="figStruct"></div>
    <table><thead><tr><th></th><th>Feature</th><th>Span (0-based)</th><th>Len</th><th>Metric</th><th>Method</th></tr></thead>
    <tbody>${rows}</tbody></table>
    <div class="flexbtns" style="margin-top:8px"><button class="btn sm" onclick="copyAllFasta('struct')">⧉ all features → FASTA</button>${csvBtns('struct')}
      <span class="lbl">click any row for its sequence / FASTA / primer</span></div>
    <div class="orfbar"><span class="lbl">ORFs · <b>${rec.orfs.length}</b> ≥40 aa</span>
      <div class="filterbox"><svg viewBox="0 0 24 24"><circle cx="11" cy="11" r="7"/><path d="m21 21-4-4"/></svg><input id="orfFilter" aria-label="Filter ORFs by strand, frame, or minimum amino acids" placeholder="filter: + / − / frame / min-aa"></div>
      <button class="btn sm" onclick="copyAllFasta('orf')">⧉ all ORFs → FASTA</button>${csvBtns('orf')}
      <span class="orfpager"><button class="btn sm" id="orfPrev">◀ prev</button>
        <span class="lbl" id="orfPageTxt"></span>
        <button class="btn sm" id="orfNext">next ▶</button></span></div>
    <table><thead><tr><th></th><th>ORF</th><th>Span (0-based)</th><th>Len</th><th>aa</th><th>Method</th></tr></thead>
    <tbody id="orfBody"></tbody></table>`
    + geneModelHtml(gm, rec.composition.length);
  state.lastRec = rec;
  glossify($("#structBody"));
  delete _gvState["figStruct"]; mountGenomeViewer($("#figStruct"), gvTracksFromRec(rec), { base:"TEagle_structure" });
  if(gm){ delete _gvState["figGene"]; mountGenomeViewer($("#figGene"), gvTracksFromGene(gm, rec.composition.length), { base:"TEagle_gene_structure" }); }
  expandCard($("#structBody").closest(".card"));
  const oi = $("#orfFilter"); if(oi) oi.oninput = e => { state.orfFilter = e.target.value; state.orfPage = 0; renderOrfPage(); };
  $("#orfPrev").onclick = () => { if(state.orfPage>0){ state.orfPage--; renderOrfPage(); } };
  $("#orfNext").onclick = () => { if((state.orfPage+1)*ORF_PP < orfFiltered().length){ state.orfPage++; renderOrfPage(); } };
  renderOrfPage();
}
/* gene structure (exon/intron/CDS) ingested from the fetched accession's annotation */
function geneModelHtml(gm, len){
  if(!gm) return "";
  const acc = (state.source && state.source.accession) || "";
  return `<div class="genebox">
    <div class="orient" style="margin-top:14px">Exons (kept in the mRNA), introns (spliced out), and the protein-coding region (<b>CDS</b>), exactly as annotated in the fetched record — displayed, not re-derived.</div>
    <div class="lbl" style="margin:12px 0 7px;font-size:12px">Gene structure — from NCBI annotation ${srcChip("NCBI", "https://www.ncbi.nlm.nih.gov/nuccore/"+encodeURIComponent(acc))}
      <span style="color:var(--faint)">· ${gm.counts.exons} exon(s) · ${gm.counts.introns} intron(s) · ${gm.counts.cds} CDS segment(s)${gm.derived_introns?" · introns derived from splice joins":""}</span></div>
    <div id="figGene"></div>
    <table><thead><tr><th>Segment</th><th>Span (0-based)</th><th>Len</th><th>Str</th><th>Note</th></tr></thead><tbody>${geneRows(gm)}</tbody></table>
    <div class="flexbtns" style="margin-top:8px"><button class="btn sm" onclick="copyGeneFasta('exon')">⧉ exons → FASTA</button>
      <button class="btn sm" onclick="copyGeneFasta('cds')">⧉ CDS → FASTA</button>${csvBtns('gene')}
      <span class="lbl">exon–intron structure ingested from the record — not re-derived</span></div>
  </div>`;
}
function _geneTag(t){ const c = t==="exon"?GENECOL.exon:t==="CDS"?GENECOL.cds:GENECOL.intron;
  return `<span class="tag" style="border-color:${c}"><span class="d" style="background:${c}"></span>${t}</span>`; }
function geneRows(gm){
  const seg = [];
  (gm.exons||[]).forEach((e,i)=>seg.push({t:"exon", s:e.start, e:e.end, str:e.strand, note:e.note?("#"+e.note):("exon "+(i+1))}));
  (gm.introns||[]).forEach((x,i)=>seg.push({t:"intron", s:x.start, e:x.end, str:x.strand, note:"intron "+(i+1)}));
  (gm.cds||[]).forEach((c)=>seg.push({t:"CDS", s:c.start, e:c.end, str:c.strand, note:c.note||"coding"}));
  seg.sort((a,b)=>a.s-b.s || (a.t>b.t?1:-1));
  return seg.map(x=>`<tr class="clickable" tabindex="0" role="button" onclick="featMenu(event, ${x.s}, ${x.e}, '${x.str}', '${x.t}_${x.s}-${x.e}')">
    <td>${_geneTag(x.t)}</td><td class="coord">${x.s}–${x.e}</td><td class="coord">${x.e-x.s}</td><td class="mono">${x.str}</td><td class="mono" style="color:var(--faint)">${x.note}</td></tr>`).join("")
    || `<tr><td colspan="5" class="mono" style="color:var(--faint);padding:10px">no exon/CDS segments in this record</td></tr>`;
}
window.copyGeneFasta = (kind) => {
  const gm = state.features; if(!gm) return;
  const arr = kind==="cds" ? (gm.cds||[]) : (gm.exons||[]);
  copyText(arr.map((x,i)=>{ const seq = x.strand==="-" ? revcomp(slice(x.start,x.end)) : slice(x.start,x.end);
    return `>${kind}_${i+1}_${x.start}-${x.end}_${x.strand||"+"}\n${seq}`; }).join("\n"), event && event.currentTarget);
};
function orfFiltered(){
  const q = (state.orfFilter||"").trim().toLowerCase();
  if(!q) return state.orfs || [];
  const num = parseInt(q, 10);
  return (state.orfs||[]).filter(o => {
    if(!isNaN(num) && String(num) === q.replace(/\D/g,"")) return o.length_aa >= num;
    const s = `${o.strand}${o.frame} ${o.strand==="+"?"+ plus":"- minus"} frame${o.frame}`.toLowerCase();
    return s.includes(q);
  });
}
function renderOrfPage(){
  const src = orfFiltered();
  const M = src.length, p = state.orfPage, a = p*ORF_PP, b = Math.min(a+ORF_PP, M);
  const N = Math.max(1, Math.ceil(M/ORF_PP));
  $("#orfBody").innerHTML = src.slice(a,b).map(o =>
    `<tr class="clickable" tabindex="0" role="button" onclick="featMenu(event, ${o.start}, ${o.end}, '${o.strand}', 'ORF_${o.strand}${o.frame}')"><td>${tagFor("ORF")}</td><td>ORF (${o.strand}${o.frame})</td>
      <td class="coord">${o.start}–${o.end}</td><td class="coord">${o.length_nt}</td>
      <td class="coord">${o.length_aa} aa</td><td class="mono" style="color:var(--faint);font-size:10.5px">6-frame ATG…stop</td></tr>`).join("")
    || `<tr><td colspan="6" class="mono" style="color:var(--faint);padding:10px">no ORFs match “${state.orfFilter||""}”</td></tr>`;
  $("#orfPageTxt").textContent = M ? `${a+1}–${b} of ${M}${state.orfFilter?" (filtered)":""} · pg ${p+1}/${N}` : "0 match";
  $("#orfPrev").disabled = p<=0;
  $("#orfNext").disabled = b>=M;
}

/* ---------- domain sequence view + copy ---------- */
window.toggleDomain = (i) => {
  const el = document.getElementById("domdet"+i); if(!el) return;
  const d = (state.domains||[])[i]; if(!d) return;
  const car = document.getElementById("dcar"+i);
  if(el.style.display === "none"){
    el.innerHTML = `
      <div class="dmeta"><a href="https://www.ebi.ac.uk/interpro/entry/pfam/${d.pfam}/" target="_blank" rel="noopener">${d.pfam} · Pfam</a> · ${d.protein.length} aa · ${d.dna.length} bp · ${d.strand} strand · nt ${d.nt[0]}–${d.nt[1]}</div>
      <div class="seqrow"><span class="slbl">DNA</span><code id="ddna${i}">${d.dna}</code></div>
      <div class="seqrow"><span class="slbl">protein</span><code id="dprot${i}">${d.protein}</code></div>
      <div class="flexbtns" style="margin-top:9px">
        <button class="btn sm" onclick="copyEl('ddna${i}',this)">copy DNA</button>
        <button class="btn sm" onclick="copyEl('dprot${i}',this)">copy protein</button>
        <button class="btn sm" onclick="copyFasta(${i},this)">copy FASTA</button>
        <button class="btn sm" onclick="designForDomain(${d.nt[0]},${d.nt[1]},'${d.domain}')">⌖ design primer</button>
      </div>`;
    el.style.display = "block"; if(car) car.textContent = "▾";
  } else { el.style.display = "none"; if(car) car.textContent = "▸"; }
};
window.copyEl = (id, btn) => copyText(document.getElementById(id)?.textContent || "", btn);
window.copyFasta = (i, btn) => { const d = state.domains[i]; copyText(`>${d.domain}_${d.pfam}_nt${d.nt[0]}-${d.nt[1]}\n${d.dna}`, btn); };
window.copyText = async (t, btn) => {
  try { await navigator.clipboard.writeText(t); if(btn){ const o = btn.textContent; btn.textContent = "copied ✓"; setTimeout(()=>btn.textContent=o, 1200); } }
  catch(e){ if(btn) btn.textContent = "copy blocked"; }
};
/* (legacy canvas map/gel removed — figures are now SVG via svgStructure/svgGel) */

/* ---------- primer design ---------- */
const PRIMER_PRESETS = {
  standard:   { pMin:150, pMax:500,  pTm:60, pMinS:18, pMaxS:27, pNum:5,  pOptS:20, pTmMin:57, pTmMax:63, pGcMin:40, pGcMax:60, pPolyX:4, pGcClamp:0 },
  qpcr:       { pMin:70,  pMax:150,  pTm:60, pMinS:18, pMaxS:24, pNum:5,  pOptS:20, pTmMin:58, pTmMax:62, pGcMin:40, pGcMax:60, pPolyX:4, pGcClamp:1 },
  highspec:   { pMin:150, pMax:500,  pTm:62, pMinS:20, pMaxS:26, pNum:8,  pOptS:22, pTmMin:60, pTmMax:64, pGcMin:45, pGcMax:60, pPolyX:3, pGcClamp:2 },
  permissive: { pMin:100, pMax:1000, pTm:58, pMinS:17, pMaxS:30, pNum:10, pOptS:20, pTmMin:52, pTmMax:65, pGcMin:30, pGcMax:70, pPolyX:5, pGcClamp:0 },
};
function applyPreset(name){
  const p = PRIMER_PRESETS[name]; if(!p) return;
  for(const k in p){ const el = $("#"+k); if(el) el.value = p[k]; }
}
function readPrimerParams(){
  const v = (id, d) => { const el = $("#"+id); return el && el.value !== "" ? +el.value : d; };
  const p = { prod_min:v("pMin",150), prod_max:v("pMax",500), opt_tm:v("pTm",60),
    min_size:v("pMinS",18), max_size:v("pMaxS",27), num_return:v("pNum",5),
    opt_size:v("pOptS",20), min_tm:v("pTmMin",57), max_tm:v("pTmMax",63),
    min_gc:v("pGcMin",40), max_gc:v("pGcMax",60), max_poly_x:v("pPolyX",4) };
  const clamp = v("pGcClamp",0); if(clamp > 0) p.gc_clamp = clamp;
  return p;
}
if($("#pPreset")) $("#pPreset").onchange = e => applyPreset(e.target.value);
if($("#pReset")) $("#pReset").onclick = () => { const s=$("#pPreset"); if(s) s.value="standard"; applyPreset("standard"); };
function resetPrimers(){                              // a failed (re)design must not leave stale, stageable pairs behind
  state.candidates = [];
  if($("#pcrStageAll")) $("#pcrStageAll").disabled = true;
  if($("#primBody")) $("#primBody").innerHTML = "";
  if($("#primMeta")) $("#primMeta").textContent = "";
}
$("#design").onclick = async () => {
  clearErr();
  if(_staleBlock()) return;
  const params = readPrimerParams();
  const btn = $("#design"); const orig = btn.textContent;
  btn.disabled = true; btn.classList.add("pending"); btn.textContent = "◴ designing…";
  try {
    const d = await api("/api/primers", { sequence: state.seq, params });
    renderPrimers(d);
    renderProvenance(d.provenance);
  } catch(e){ resetPrimers(); err("Primer design failed — " + e.message); }
  finally { btn.disabled = false; btn.classList.remove("pending"); btn.textContent = orig; }
};
window.designForDomain = async (start, end, label) => {
  clearErr();
  if(_staleBlock()) return;                          // same staleness guard as #design / #runpcr
  const params = readPrimerParams();
  scrollToCard($("#design").closest(".card"));
  $("#primMeta").textContent = `domain ${label} · designing…`;
  try {
    const d = await api("/api/primers", { sequence: state.seq, params, included: [start, Math.max(60, end - start)] });
    renderPrimers(d, `Domain-specific PCR — primers placed within the ${label} domain (nt ${start}–${end}). A conserved-domain primer may amplify multiple TE families or copies; confirm with in-silico PCR before use.`);
    renderProvenance(d.provenance);
    $("#primMeta").textContent = `${d.candidates.length} pairs · domain ${label}`;
  } catch(e){ resetPrimers(); err("Domain primer design failed — " + e.message); }
};

function hl3(seq){ // underline the true 3' terminal 5 bases (both primers are written 5'->3')
  const n = 5;
  const tip = "3′ terminus — these bases must match the template exactly for the primer to extend; they govern PCR specificity";
  return seq.slice(0, -n) + `<span class="p3" title="${tip}">${seq.slice(-n)}</span>`;
}
function renderPrimers(d, warn){
  expandCard($("#primBody").closest(".card"));
  state.candidates = d.candidates || [];
  if($("#pcrStageAll")) $("#pcrStageAll").disabled = !state.candidates.length;
  const warnHtml = warn ? `<div class="pcrwarn">⚠ ${warn}</div>` : "";
  if(!d.candidates.length){ $("#primBody").innerHTML = warnHtml + `<div class="empty">Primer3 returned no pairs — ${d.explain_pair||"loosen constraints"}.</div>`; return; }
  $("#primBody").innerHTML = warnHtml + `<div class="flexbtns" style="margin-bottom:10px"><span class="lbl">${d.candidates.length} primer pair(s)</span>${csvBtns('primers')}</div>` + d.candidates.map((c,i)=>`
    <div class="oligo">
      <div class="oh"><span class="id">${c.id}</span>
        <span class="tag t-good" title="Primer3's overall penalty for this pair — lower is better; it rises as the primers depart from your target Tm, size, and GC"><span class="d" style="background:var(--good)"></span>penalty ${c.penalty}</span>
        <span class="meta" style="margin-left:auto;font-family:var(--mono);font-size:10.5px;color:var(--faint)">product ${c.product_size} bp</span></div>
      <div class="ob">
        <div class="strand"><span class="s">FWD</span><span class="seq">${hl3(c.left_seq)}</span><span class="m" title="Melting temperature (°C) — the temperature at which half the primer–template duplex separates; the two primers are matched so both anneal together">Tm ${c.left_tm}</span><span class="m" title="Percent G+C content of the primer; ~40–60% is typical">GC ${c.left_gc}%</span></div>
        <div class="strand"><span class="s">REV</span><span class="seq">${hl3(c.right_seq)}</span><span class="m" title="Melting temperature (°C) — matched to the forward primer">Tm ${c.right_tm}</span><span class="m" title="Percent G+C content of the primer; ~40–60% is typical">GC ${c.right_gc}%</span></div>
        <div class="flexbtns" style="margin-top:9px">
          <button class="btn sm" onclick='selectPair(${JSON.stringify(c).replace(/'/g,"&#39;")})'>→ send to in-silico PCR</button>
          <button class="btn sm" onclick="copyText('>${c.id}_FWD\\n${c.left_seq}\\n>${c.id}_REV\\n${c.right_seq}', this)">⧉ copy pair</button>
          <span class="lbl" title="The 3′ end is where the polymerase extends; TEagle's in-silico PCR requires an exact match here, so it is underlined.">3′ terminal 5 bases underlined — the specificity-determining end</span></div>
      </div>
    </div>`).join("");
  $("#primMeta").textContent = `${d.candidates.length} pairs · Primer3`;
}
/* ---------- in-silico PCR engine: staged primer-pair manager (see / remove / reorder) ---------- */
function pcrKey(c){ return c.left_seq + "|" + c.right_seq; }
function renderPcrQueue(){
  const host = $("#pcrQueue"); if(!host) return;
  const n = state.pcrPairs.length;
  if($("#pcrCount")) $("#pcrCount").textContent = `${n} loaded`;
  if($("#runpcr")) $("#runpcr").disabled = !n;
  if($("#pcrClear")) $("#pcrClear").disabled = !n;
  if(!n){ host.innerHTML = `<div class="empty sm">No pairs loaded. Design primers (panel 04), then “→ send to in-silico PCR”, or stage all below.</div>`; return; }
  host.innerHTML = state.pcrPairs.map((c,i)=>`
    <div class="pcrrow" draggable="true" data-i="${i}"
         ondragstart="pcrDragStart(event,${i})" ondragover="pcrDragOver(event,${i})" ondragleave="pcrDragLeave(event)" ondrop="pcrDrop(event,${i})" ondragend="pcrDragLeave(event)">
      <span class="grip" title="drag to reorder" aria-hidden="true">⠿</span>
      <span class="ln" title="gel lane ${i+1}">L${i+1}</span>
      <span class="pid">${c.id||("pair"+(i+1))}</span>
      <span class="pseq mono" title="forward primer 5′→3′">F ${c.left_seq}</span>
      <span class="pseq mono" title="reverse primer 5′→3′">R ${c.right_seq}</span>
      <span class="psz">${c.product_size||"?"} bp</span>
      <span class="pcrctl">
        <button class="btn xs" title="move up (earlier lane)" ${i===0?"disabled":""} onclick="movePcrPair(${i},-1)" aria-label="move pair up">↑</button>
        <button class="btn xs" title="move down (later lane)" ${i===n-1?"disabled":""} onclick="movePcrPair(${i},1)" aria-label="move pair down">↓</button>
        <button class="btn xs" title="remove this pair" onclick="removePcrPair(${i})" aria-label="remove pair">✕</button>
      </span>
    </div>`).join("");
}
let _pcrDrag = null;
window.pcrDragStart = (e,i) => { _pcrDrag = i; if(e.dataTransfer){ e.dataTransfer.effectAllowed = "move"; try{ e.dataTransfer.setData("text/plain", String(i)); }catch(_){} } };
window.pcrDragOver = (e,i) => { e.preventDefault(); if(e.dataTransfer) e.dataTransfer.dropEffect = "move"; const r=e.currentTarget; if(r) r.classList.add("dragover"); };
window.pcrDragLeave = (e) => { const r=e.currentTarget; if(r) r.classList.remove("dragover"); };
window.pcrDrop = (e,i) => { e.preventDefault(); const from = _pcrDrag; _pcrDrag = null;
  if(from == null || from === i) return renderPcrQueue();
  const a = state.pcrPairs; const [m] = a.splice(from,1); a.splice(i,0,m); renderPcrQueue(); };
window.addPcrPair = (c) => {
  const k = pcrKey(c);
  if(!state.pcrPairs.some(p => pcrKey(p) === k)) state.pcrPairs.push(c);
  renderPcrQueue();
  $("#pcrHint").textContent = `${state.pcrPairs.length} pair(s) loaded · run to search`;
  expandCard($("#runpcr").closest(".card"));
};
window.removePcrPair = (i) => { state.pcrPairs.splice(i,1); renderPcrQueue(); };
window.movePcrPair = (i,d) => { const j=i+d, a=state.pcrPairs; if(j<0||j>=a.length) return; [a[i],a[j]]=[a[j],a[i]]; renderPcrQueue(); };
window.selectPair = (c) => addPcrPair(c);                     // "→ send to in-silico PCR" stages the pair
if($("#pcrClear")) $("#pcrClear").onclick = () => { state.pcrPairs = []; renderPcrQueue(); };
if($("#pcrStageAll")) $("#pcrStageAll").onclick = () => {
  (state.candidates||[]).forEach(c => { const k=pcrKey(c); if(!state.pcrPairs.some(p=>pcrKey(p)===k)) state.pcrPairs.push(c); });
  renderPcrQueue(); expandCard($("#runpcr").closest(".card"));
};

/* ---------- in-silico PCR: run the staged pairs (order = lane order) ---------- */
$("#runpcr").onclick = async () => {
  clearErr();
  if(_staleBlock()) return;
  const pairs = state.pcrPairs || [];
  if(!pairs.length){ err("Load at least one primer pair into the engine first (design primers, then “→ send to in-silico PCR”)."); return; }
  const btn = $("#runpcr"), o = btn.textContent; btn.disabled = true; btn.classList.add("pending"); btn.textContent = "◴ running…";
  expandCard($("#runpcr").closest(".card"));
  const pv = (id, d) => { const el = $("#"+id); const v = (el && el.value !== "") ? +el.value : NaN; return Number.isFinite(v) ? v : d; };  // a blank field falls back to the default (not 0)
  const bg = $("#bg").value, p = { max_mm: pv("pcrMM",2), tp: pv("pcrTP",5), prod_min: pv("pcrPmin",70), prod_max: pv("pcrPmax",1000) };
  const runOne = c => api("/api/pcr", { sequence: state.seq, background: bg, fwd: c.left_seq, rev: c.right_seq,
    target_span: [c.left_pos[0], c.right_pos[1]], params: p });
  try {
    if(pairs.length === 1){                                   // single pair → detailed single-lane view
      const d = await runOne(pairs[0]);
      renderPCR(d); renderProvenance(d.provenance);
    } else {                                                  // multiple → multi-lane gel, one lane per staged pair
      $("#pcrBody").innerHTML = `<div class="mini">running in-silico PCR for ${pairs.length} loaded primer pairs…</div>`;
      const lanes = [], amps = [], provs = [];
      pairs.forEach((c, i) => { c.__lane = `P${i + 1}`; });     // unique lane names in staged order (avoid duplicate ids)
      for(let i = 0; i < pairs.length; i++){
        const c = pairs[i], lane = `P${i + 1}`;
        try {
          const d = await runOne(c);
          lanes.push({ label: lane, amplicons: d.amplicons });
          d.amplicons.forEach(a => amps.push({ ...a, pair: lane }));
          if(d.provenance) provs.push(d.provenance);           // every pair is its own sealed manifest — keep them all
        } catch(e){ lanes.push({ label: lane + " ✕", amplicons: [], error: true }); }
      }
      state.lastPcr = { lanes, amplicons: amps, provenances: provs };
      renderPCRMulti(lanes, amps);
      if(provs.length){
        renderProvenance(provs[0]);
        if(provs.length > 1){                                  // disclose that each pair was sealed independently (panel shows P1)
          const note = document.createElement("div");
          note.className = "mini"; note.style.marginBottom = "8px";
          note.textContent = `${provs.length} primer pairs were run — each sealed as its own manifest (showing P1; all ${provs.length} recorded this session).`;
          $("#provBody").prepend(note);
        }
      }
    }
  } catch(e){ err("In-silico PCR failed — " + e.message); }
  finally { btn.disabled = false; btn.classList.remove("pending"); btn.textContent = o; }
};
function renderPCRMulti(lanes, amps){
  expandCard($("#pcrBody").closest(".card"));
  const onN = amps.filter(a=>a.on_target).length;
  const rows = amps.map(a => `<tr class="clickable" tabindex="0" role="button" onclick="featMenu(event, ${a.start}, ${a.end}, '+', 'amplicon_${a.pair}')">
    <td class="mono">${a.pair}</td><td class="mono" style="font-size:11px">${esc(a.source)}</td>
    <td class="coord">${a.start.toLocaleString()}–${a.end.toLocaleString()}</td><td class="coord">${a.length}</td><td class="coord">${a.fwd_mm}/${a.rev_mm}</td>
    <td>${a.on_target?'<span class="tag t-good"><span class="d" style="background:var(--good)"></span>on-target</span>':'<span class="tag t-bad"><span class="d" style="background:var(--bad)"></span>off-target</span>'}</td></tr>`).join("");
  $("#pcrBody").innerHTML = `
    <div id="figGel"></div>
    ${amps.length?`<table><thead><tr><th>Pair</th><th>Source</th><th>Coords</th><th>Len</th><th>Mism F/R</th><th>Call</th></tr></thead><tbody>${rows}</tbody></table>
      <div class="flexbtns" style="margin-top:8px"><button class="btn sm" onclick="copyAmplicons()">⧉ amplicons → FASTA</button>${csvBtns('pcr')}</div>`
      :`<div class="empty">No amplicon predicted for any pair under the criteria.</div>`}
    <div class="mini" style="margin-top:10px">${lanes.length} lanes · one per primer pair · ${onN} on-target band(s) · ladder lane “L”. Not a claim of experimental specificity.</div>`;
  glossify($("#pcrBody"));
  mountFigure($("#figGel"), b => svgGel({ lanes }, b), "TEagle_gel_multilane", { cls:"gel", modes:["dark","white","uv","mono"], viewH:470 });
}
function renderPCR(d){
  const on = d.amplicons.filter(a=>a.on_target).length, off = d.amplicons.length-on;
  const rows = d.amplicons.map((a,i)=>`
    <tr><td class="coord">${i+1}</td><td class="mono" style="font-size:11px">${esc(a.source)}</td>
      <td class="coord">${a.start.toLocaleString()}–${a.end.toLocaleString()}</td>
      <td class="coord">${a.length}</td><td class="coord">${a.fwd_mm}/${a.rev_mm}</td>
      <td>${a.on_target?'<span class="tag t-good"><span class="d" style="background:var(--good)"></span>on-target</span>':'<span class="tag t-bad"><span class="d" style="background:var(--bad)"></span>off-target</span>'}</td></tr>`).join("");
  const bgs = d.backgroundsSearched.map(b=>`<span class="tag t-good"><span class="d" style="background:var(--good)"></span>${b}</span>`).join(" ");
  const nots = d.notSearched.map(b=>`<div class="notrun"><span class="x">✕</span>${b} — not run</div>`).join("");
  const verdict = off===0
    ? `No additional amplicon was predicted in the searched sequences under the specified mismatch, orientation and product-length criteria. This is not a claim of experimental specificity.`
    : `${off} additional (off-target) amplicon(s) predicted in the searched sequences under the current criteria.`;
  $("#pcrBody").innerHTML = `
    <div id="figGel"></div>
    <div style="margin:6px 0 12px">${bgs}</div>
    ${d.amplicons.length?`<table><thead><tr><th>#</th><th>Source</th><th>Coords</th><th>Len</th><th>Mism F/R</th><th>Call</th></tr></thead><tbody>${rows}</tbody></table>
      <div class="flexbtns" style="margin-top:8px"><button class="btn sm" onclick="copyAmplicons()">⧉ amplicons → FASTA</button>${csvBtns('pcr')}</div>`
      :`<div class="empty">No amplicon predicted under the criteria.</div>`}
    <div class="explain" style="margin-top:13px">${verdict}</div>
    <div style="margin-top:11px">${nots}</div>
    <div class="mini">criteria · ≤${d.criteria.max_mismatch} mismatches · strict 3′ (last ${d.criteria.three_prime_strict}) · product ${d.criteria.product_size[0]}–${d.criteria.product_size[1]} bp</div>`;
  state.lastPcr = d;
  glossify($("#pcrBody"));
  mountFigure($("#figGel"), b => svgGel(d, b), "TEagle_gel", { cls:"gel", modes:["dark","white","uv","mono"], viewH:470 });
  expandCard($("#pcrBody").closest(".card"));
}

/* ---------- provenance ---------- */
function renderProvenance(m){
  expandCard($("#provBody").closest(".card"));
  const sw = m.software.map(s=>`<div class="li">${esc(s.name)}<span>${esc(s.version)}</span></div>`).join("");
  const pr = Object.entries(m.parameters||{}).map(([k,v])=>`<div class="li">${esc(k)}<span class="${String(v).length>18?'hash':''}">${esc(v==null?"—":String(v))}</span></div>`).join("") || `<div class="li">defaults<span>—</span></div>`;  // the manifest is the reproducibility record — show EVERY sealed parameter
  const db = (m.databases||[]).map(d=>`<div class="li">${esc(d.name)}<span class="hash">${esc(d.version || (d.sha256?d.sha256.slice(0,12)+"…":(d.file||"—")))}</span></div>`).join("");
  const dbBlock = db ? `<div class="mblock"><div class="t">Databases</div>${db}</div>` : "";
  const nr = m.notRun.map(n=>`<div class="notrun"><span class="x">✕</span>${esc(n)}</div>`).join("");
  const rf = (m.references||[]).map(r=>`<div class="refitem"><b>${r.name}</b> — ${r.citation} ${r.doi?`<a href="https://doi.org/${r.doi}" target="_blank" rel="noopener">doi:${r.doi}</a>`:""}${r.license?` · ${r.license}`:""}</div>`).join("");
  const refBlock = rf ? `<div class="refblock"><div class="lbl" style="margin-bottom:8px">Database &amp; tool references (source-verified)</div>${rf}</div>` : "";
  $("#provBody").innerHTML = `
    <div class="manifest">
      <div class="mblock"><div class="t">Input</div>
        <div class="li">id<span>${esc(m.input.id)}</span></div>
        <div class="li">length<span>${m.input.length.toLocaleString()} bp</span></div>
        <div class="li">sha256<span class="hash">${m.input.sha256.slice(0,16)}…</span></div>
        <div class="li">run type<span>${esc(m.runType)}</span></div></div>
      <div class="mblock"><div class="t">Software</div>${sw}</div>
      <div class="mblock"><div class="t">Parameters</div>${pr}</div>
      ${dbBlock}
      <div class="mblock"><div class="t">Environment</div>
        <div class="li">os<span style="font-size:9.5px">${(m.environment.os||"").slice(0,22)}</span></div>
        <div class="li">python<span>${m.environment.python}</span></div>
        <div class="li">created<span style="font-size:9.5px">${(m.createdUtc||"").replace("T"," ").replace("+00:00","Z")}</span></div>
        <div class="li">manifest<span class="hash">${m.manifestSha256.slice(0,14)}…</span></div></div>
    </div>
    <div style="margin-top:12px">${nr}</div>${refBlock}`;
}
/* SVG figures are responsive (width:100%); no manual resize redraw needed */
