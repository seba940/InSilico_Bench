"""Auto Primer Design tab using primer3.

Modes
-----
1. CDS Up & Down  – amplify entire CDS from flanking regions
2. CDS Internal   – amplify portion inside CDS (genotyping / cloning)
3. RT-PCR / qPCR  – short amplicons, tight Tm/GC constraints
4. Junction Check – one external primer + one internal primer for integration check
"""
import tkinter as tk
from tkinter import ttk, messagebox

from core.models import Primer
from gui.gui_components import UIHelper, SearchableCombobox, SnapGeneViewer

try:
    import primer3
    _P3_OK = True
except ImportError:
    _P3_OK = False


def _gc(seq):
    if not seq:
        return 0.0
    s = seq.upper()
    return round(100.0 * (s.count("G") + s.count("C")) / len(s), 1)


def _calc_tm(seq):
    try:
        try:
            from Bio.SeqUtils.MeltingTemp import Tm_NN
        except ImportError:
            from Bio.SeqUtils import Tm_NN
        return round(Tm_NN(seq, Na=50, dnac1=250, dnac2=250), 1)
    except Exception:
        return 0.0


_COMMON_GLOBAL = {
    "PRIMER_OPT_SIZE": 20,
    "PRIMER_MIN_SIZE": 18,
    "PRIMER_MAX_SIZE": 27,
    "PRIMER_OPT_TM": 60.0,
    "PRIMER_MIN_TM": 57.0,
    "PRIMER_MAX_TM": 63.0,
    "PRIMER_MIN_GC": 40.0,
    "PRIMER_MAX_GC": 60.0,
    "PRIMER_MAX_POLY_X": 4,
    "PRIMER_SALT_MONOVALENT": 50.0,
    "PRIMER_DNA_CONC": 250.0,
    "PRIMER_NUM_RETURN": 5,
    "PRIMER_EXPLAIN_FLAG": 1,
}


class PrimerDesignTab:
    """primer3-based automatic primer designer."""

    def __init__(self, parent, lib, refresh_callback):
        self.parent = parent
        self.lib = lib
        self.refresh_cb = refresh_callback
        self._results = []   # list of dicts from last run, one per pair
        self._template_seq = ""
        self._local_start = 0   # gene start relative to viewer subseq
        self._local_end = 0     # gene end relative to viewer subseq
        self._flank = 200
        self.setup_ui()

    # ── UI ───────────────────────────────────────────────────────────────

    def setup_ui(self):
        if not _P3_OK:
            ttk.Label(self.parent, text="primer3-py가 설치되지 않았습니다.\npip install primer3-py",
                      foreground="red", font=("TkDefaultFont", 11)).pack(expand=True)
            return

        outer = ttk.Frame(self.parent)
        outer.pack(fill="both", expand=True)

        # Left panel (controls)
        left = ttk.Frame(outer, width=400)
        left.pack(side="left", fill="y", padx=4, pady=4)
        left.pack_propagate(False)

        # Right panel (viewer)
        right = ttk.Frame(outer)
        right.pack(side="left", fill="both", expand=True, padx=4, pady=4)

        self._build_left(left)
        self._build_right(right)

    def _build_left(self, f):
        # ── 1. Target ────────────────────────────────────────────────────
        src = ttk.LabelFrame(f, text="1. Target", padding=6)
        src.pack(fill="x", pady=3)

        r0 = ttk.Frame(src); r0.pack(fill="x", pady=2)
        ttk.Label(r0, text="Template:", width=16).pack(side="left")
        self.tmp_c = SearchableCombobox(r0)
        self.tmp_c.pack(side="left", fill="x", expand=True, padx=2)
        self.tmp_c.bind("<<ComboboxSelected>>", lambda e: self._on_template_change())

        r1 = ttk.Frame(src); r1.pack(fill="x", pady=2)
        ttk.Label(r1, text="Feature / Gene:", width=16).pack(side="left")
        self.feat_c = ttk.Combobox(r1, state="readonly")
        self.feat_c.pack(side="left", fill="x", expand=True, padx=2)
        self.feat_c.bind("<<ComboboxSelected>>", lambda e: self._on_feature_change())

        r2 = ttk.Frame(src); r2.pack(fill="x", pady=2)
        ttk.Label(r2, text="Manual Start:", width=16).pack(side="left")
        self.coord_s = ttk.Entry(r2, width=10); self.coord_s.pack(side="left", padx=2)
        ttk.Label(r2, text="End:").pack(side="left")
        self.coord_e = ttk.Entry(r2, width=10); self.coord_e.pack(side="left", padx=2)

        r2b = ttk.Frame(src); r2b.pack(fill="x", pady=2)
        ttk.Label(r2b, text="Strand:", width=16).pack(side="left")
        self.strand_v = tk.StringVar(value="+")
        ttk.Radiobutton(r2b, text="+", variable=self.strand_v, value="+").pack(side="left")
        ttk.Radiobutton(r2b, text="-", variable=self.strand_v, value="-").pack(side="left")

        self.target_info = ttk.Label(src, text="", foreground="gray",
                                     font=("TkDefaultFont", 8, "italic"))
        self.target_info.pack(anchor="w")

        # ── 2. Mode ──────────────────────────────────────────────────────
        mode_f = ttk.LabelFrame(f, text="2. Design Mode", padding=6)
        mode_f.pack(fill="x", pady=3)

        self.mode_v = tk.StringVar(value="updown")
        modes = [
            ("CDS Up & Down (전체 증폭)",            "updown"),
            ("CDS Internal (내부 증폭 / Genotyping)", "internal"),
            ("RT-PCR / qPCR (정량 분석)",             "rtpcr"),
            ("Junction Check 5' (HR 5' 검증)",       "junc5"),
            ("Junction Check 3' (HR 3' 검증)",       "junc3"),
        ]
        for text, val in modes:
            ttk.Radiobutton(mode_f, text=text, variable=self.mode_v, value=val).pack(anchor="w")

        # ── 3. Parameters ────────────────────────────────────────────────
        par_f = ttk.LabelFrame(f, text="3. Parameters", padding=6)
        par_f.pack(fill="x", pady=3)

        rp0 = ttk.Frame(par_f); rp0.pack(fill="x", pady=1)
        ttk.Label(rp0, text="Flanking (bp):", width=20).pack(side="left")
        self.flank_sp = ttk.Spinbox(rp0, from_=50, to=1000, width=7)
        self.flank_sp.set(200); self.flank_sp.pack(side="left", padx=2)

        rp1 = ttk.Frame(par_f); rp1.pack(fill="x", pady=1)
        ttk.Label(rp1, text="Internal size (bp):", width=20).pack(side="left")
        self.internal_sp = ttk.Spinbox(rp1, from_=30, to=200, width=7)
        self.internal_sp.set(80); self.internal_sp.pack(side="left", padx=2)
        ttk.Label(rp1, text="(Junction Check only)", foreground="gray",
                  font=("TkDefaultFont", 7)).pack(side="left")

        rp2 = ttk.Frame(par_f); rp2.pack(fill="x", pady=1)
        ttk.Label(rp2, text="Opt Tm (°C):", width=20).pack(side="left")
        self.tm_opt = ttk.Spinbox(rp2, from_=50, to=75, width=7, increment=0.5)
        self.tm_opt.set(60.0); self.tm_opt.pack(side="left", padx=2)

        rp3 = ttk.Frame(par_f); rp3.pack(fill="x", pady=1)
        ttk.Label(rp3, text="Min/Max Tm:", width=20).pack(side="left")
        self.tm_min = ttk.Spinbox(rp3, from_=45, to=70, width=6, increment=0.5)
        self.tm_min.set(57.0); self.tm_min.pack(side="left", padx=2)
        ttk.Label(rp3, text="–").pack(side="left")
        self.tm_max = ttk.Spinbox(rp3, from_=55, to=80, width=6, increment=0.5)
        self.tm_max.set(63.0); self.tm_max.pack(side="left", padx=2)

        # ── Run button ───────────────────────────────────────────────────
        ttk.Button(f, text="⚡ Design Primers", command=self.run).pack(pady=6, fill="x")

        # ── 4. Results table ────────────────────────────────────────────
        res_f = ttk.LabelFrame(f, text="4. Results (Top 5)", padding=6)
        res_f.pack(fill="both", expand=True, pady=3)

        cols = ("rank", "fwd", "rev", "tm_f", "tm_r", "gc_f", "gc_r", "size", "penalty")
        heads = ("#", "Fwd Seq", "Rev Seq", "Tm_F", "Tm_R", "GC%_F", "GC%_R", "Size", "Penalty")
        fr_t, self.tree = UIHelper.create_scrolled_tree(res_f, cols, heads)
        fr_t.pack(fill="both", expand=True)
        widths = [25, 160, 160, 52, 52, 52, 52, 52, 60]
        for col, w in zip(cols, widths):
            self.tree.column(col, width=w)
        self.tree.bind("<<TreeviewSelect>>", self._on_result_select)

        # Save controls
        sav_f = ttk.Frame(res_f); sav_f.pack(fill="x", pady=3)
        ttk.Label(sav_f, text="Name prefix:").pack(side="left")
        self.name_pfx = ttk.Entry(sav_f, width=16)
        self.name_pfx.pack(side="left", padx=3)
        ttk.Button(sav_f, text="Save Selected", command=self.save_selected).pack(side="left", padx=2)
        ttk.Button(sav_f, text="Save All", command=self.save_all).pack(side="left")

    def _build_right(self, f):
        ttk.Label(f, text="Sequence Preview", font=("TkDefaultFont", 9, "bold")).pack(anchor="w")
        self.viewer = SnapGeneViewer(f, lib_manager=self.lib)
        self.viewer.pack(fill="both", expand=True)

    # ── Source refresh ────────────────────────────────────────────────────

    def refresh_sources(self):
        if not _P3_OK:
            return
        all_tmpl = []
        for k, v in self.lib.templates.items():
            all_tmpl.append(f"[{v.topology}] {k}")
        for k, v in self.lib.digests.items():
            all_tmpl.append(f"[Digest:{v.topology}] {k}")
        for k, v in self.lib.recombinants.items():
            all_tmpl.append(f"[HR:{v.topology}] {k}")
        for k, v in self.lib.amplicons.items():
            all_tmpl.append(f"[Amplicon] {k}")
        self.tmp_c['values'] = sorted(all_tmpl)

    # ── Event handlers ────────────────────────────────────────────────────

    def _resolve_template(self, label):
        """Return (name, SequenceItem) for a display label."""
        if label.startswith("[Digest:"):
            name = label.split("] ", 1)[1]
            return name, self.lib.digests.get(name)
        if label.startswith("[HR:"):
            name = label.split("] ", 1)[1]
            return name, self.lib.recombinants.get(name)
        if label.startswith("[Amplicon]"):
            name = label.split("] ", 1)[1]
            return name, self.lib.amplicons.get(name)
        # [topology] name
        name = label.split("] ", 1)[-1] if "] " in label else label
        return name, self.lib.templates.get(name)

    def _on_template_change(self):
        label = self.tmp_c.get()
        _, item = self._resolve_template(label)
        if item is None:
            return
        self._template_seq = item.sequence
        feats = getattr(item, "features", [])
        gene_names = sorted({f.get("label", "") for f in feats if f.get("label")})
        self.feat_c['values'] = gene_names[:500]
        self.feat_c.set("")
        self.target_info.config(text=f"{len(item.sequence):,} bp ({item.topology})",
                                foreground="gray")
        self._render_viewer([], [])

    def _on_feature_change(self):
        label = self.tmp_c.get()
        _, item = self._resolve_template(label)
        if item is None:
            return
        feat_name = self.feat_c.get()
        for feat in getattr(item, "features", []):
            if feat.get("label") == feat_name:
                s, e = int(feat.get("start", 0)), int(feat.get("end", 0))
                strand = "+" if feat.get("strand", 1) >= 0 else "-"
                self.coord_s.delete(0, tk.END); self.coord_s.insert(0, str(s))
                self.coord_e.delete(0, tk.END); self.coord_e.insert(0, str(e))
                self.strand_v.set(strand)
                self.target_info.config(
                    text=f"{feat_name}  {s}–{e} ({strand})  {e-s} bp",
                    foreground="darkgreen")
                break

    # ── Core design ───────────────────────────────────────────────────────

    def _get_coords(self):
        try:
            gs = int(self.coord_s.get())
            ge = int(self.coord_e.get())
            flank = int(self.flank_sp.get())
            return gs, ge, flank
        except ValueError:
            raise ValueError("Start/End/Flanking에 올바른 숫자를 입력하세요.")

    def _build_args(self, gene_seq, gene_len, local_s, local_e, mode, flank):
        """Return (seq_args, global_args) for primer3 based on mode."""
        try:
            tm_opt = float(self.tm_opt.get())
            tm_min = float(self.tm_min.get())
            tm_max = float(self.tm_max.get())
            internal_size = int(self.internal_sp.get())
        except ValueError:
            tm_opt, tm_min, tm_max = 60.0, 57.0, 63.0
            internal_size = 80

        g = dict(_COMMON_GLOBAL)
        g["PRIMER_OPT_TM"] = tm_opt
        g["PRIMER_MIN_TM"] = tm_min
        g["PRIMER_MAX_TM"] = tm_max

        if mode == "updown":
            # Amplify entire CDS: primers sit in flanking regions
            min_prod = gene_len + 60
            max_prod = gene_len + max(flank, 300)
            seq_args = {
                "SEQUENCE_ID": "target",
                "SEQUENCE_TEMPLATE": gene_seq,
                "SEQUENCE_TARGET": [[local_s, gene_len]],
            }
            g["PRIMER_PRODUCT_SIZE_RANGE"] = [[min_prod, max_prod]]

        elif mode == "internal":
            min_prod = min(200, max(100, gene_len // 4))
            max_prod = min(1000, gene_len)
            if min_prod >= max_prod:
                min_prod = max(80, max_prod - 100)
            seq_args = {
                "SEQUENCE_ID": "target",
                "SEQUENCE_TEMPLATE": gene_seq,
                "SEQUENCE_INCLUDED_REGION": [local_s, gene_len],
            }
            g["PRIMER_PRODUCT_SIZE_RANGE"] = [[min_prod, max_prod]]

        elif mode == "rtpcr":
            seq_args = {
                "SEQUENCE_ID": "target",
                "SEQUENCE_TEMPLATE": gene_seq,
                "SEQUENCE_INCLUDED_REGION": [local_s, gene_len],
            }
            g["PRIMER_PRODUCT_SIZE_RANGE"] = [[75, 250]]
            g["PRIMER_OPT_TM"] = 60.0
            g["PRIMER_MIN_TM"] = 58.0
            g["PRIMER_MAX_TM"] = 62.0
            g["PRIMER_MIN_GC"] = 40.0
            g["PRIMER_MAX_GC"] = 60.0

        elif mode == "junc5":
            # Left primer in upstream flank, right primer in first part of gene
            right_region_start = local_s
            right_region_len = min(internal_size, gene_len)
            min_prod = right_region_len + 30
            max_prod = right_region_len + flank
            seq_args = {
                "SEQUENCE_ID": "target",
                "SEQUENCE_TEMPLATE": gene_seq,
                "SEQUENCE_PRIMER_PAIR_OK_REGION_LIST": [
                    [0, local_s, right_region_start, right_region_len]
                ],
            }
            g["PRIMER_PRODUCT_SIZE_RANGE"] = [[min_prod, max_prod]]

        elif mode == "junc3":
            # Left primer in last part of gene, right primer in downstream flank
            left_region_start = max(0, local_e - internal_size)
            left_region_len = local_e - left_region_start
            min_prod = left_region_len + 30
            max_prod = left_region_len + flank
            seq_args = {
                "SEQUENCE_ID": "target",
                "SEQUENCE_TEMPLATE": gene_seq,
                "SEQUENCE_PRIMER_PAIR_OK_REGION_LIST": [
                    [left_region_start, left_region_len, local_e, len(gene_seq) - local_e]
                ],
            }
            g["PRIMER_PRODUCT_SIZE_RANGE"] = [[min_prod, max_prod]]

        else:
            raise ValueError(f"Unknown mode: {mode}")

        return seq_args, g

    def _parse_results(self, p3_out, offset):
        """Extract top-5 primer pairs from primer3 output dict."""
        results = []
        for i in range(5):
            lk = f"PRIMER_LEFT_{i}"
            rk = f"PRIMER_RIGHT_{i}"
            if lk not in p3_out:
                break
            lpos = p3_out[lk]          # [start, length]
            rpos = p3_out[rk]          # [3prime_end, length]
            fwd_seq = p3_out.get(f"PRIMER_LEFT_{i}_SEQUENCE", "")
            rev_seq = p3_out.get(f"PRIMER_RIGHT_{i}_SEQUENCE", "")
            penalty = p3_out.get(f"PRIMER_PAIR_{i}_PENALTY", 0.0)
            size = p3_out.get(f"PRIMER_PAIR_{i}_PRODUCT_SIZE", 0)
            tm_f = p3_out.get(f"PRIMER_LEFT_{i}_TM", _calc_tm(fwd_seq))
            tm_r = p3_out.get(f"PRIMER_RIGHT_{i}_TM", _calc_tm(rev_seq))
            # absolute positions on the full template
            abs_lstart = lpos[0] + offset
            abs_rend = rpos[0] + offset      # 3'-end on + strand
            abs_rstart = abs_rend - rpos[1] + 1
            results.append({
                "rank": i + 1,
                "fwd_seq": fwd_seq.upper(),
                "rev_seq": rev_seq.upper(),
                "tm_f": round(tm_f, 1),
                "tm_r": round(tm_r, 1),
                "gc_f": _gc(fwd_seq),
                "gc_r": _gc(rev_seq),
                "size": size,
                "penalty": round(penalty, 3),
                "lpos": lpos,           # [start, len] relative to gene_seq
                "rpos": rpos,           # [3p_end, len] relative to gene_seq
                "abs_lstart": abs_lstart,
                "abs_rstart": abs_rstart,
                "abs_rend": abs_rend,
            })
        return results

    def run(self):
        label = self.tmp_c.get()
        if not label:
            messagebox.showwarning("Warning", "Template을 선택하세요.")
            return
        _, item = self._resolve_template(label)
        if item is None:
            messagebox.showwarning("Warning", "Template을 찾을 수 없습니다.")
            return

        try:
            gene_start, gene_end, flank = self._get_coords()
        except ValueError as exc:
            messagebox.showerror("Error", str(exc))
            return

        full_seq = item.sequence
        gene_len = gene_end - gene_start
        if gene_len <= 0:
            messagebox.showerror("Error", "유효한 Start/End 좌표를 입력하세요.")
            return

        # Extract subseq with flanking regions
        sub_start = max(0, gene_start - flank)
        sub_end = min(len(full_seq), gene_end + flank)
        gene_seq = full_seq[sub_start:sub_end]
        local_s = gene_start - sub_start
        local_e = gene_end - sub_start
        self._local_start = local_s
        self._local_end = local_e
        self._flank = flank

        if self.strand_v.get() == "-":
            from core.utils import get_rc
            gene_seq = get_rc(gene_seq)
            gene_len_seq = len(gene_seq)
            local_s_rc = gene_len_seq - local_e
            local_e_rc = gene_len_seq - local_s
            local_s, local_e = local_s_rc, local_e_rc

        mode = self.mode_v.get()
        try:
            seq_args, global_args = self._build_args(
                gene_seq, gene_len, local_s, local_e, mode, flank)
        except ValueError as exc:
            messagebox.showerror("Error", str(exc))
            return

        try:
            p3_out = primer3.design_primers(seq_args, global_args)
        except Exception as exc:
            messagebox.showerror("primer3 Error", str(exc))
            return

        n_pairs = p3_out.get("PRIMER_PAIR_NUM_RETURNED", 0)
        if n_pairs == 0:
            exp = p3_out.get("PRIMER_PAIR_EXPLAIN", p3_out.get("PRIMER_LEFT_EXPLAIN", ""))
            messagebox.showinfo("No Results",
                                f"조건에 맞는 프라이머 쌍을 찾지 못했습니다.\n{exp}")
            self._results = []
            self.tree.delete(*self.tree.get_children())
            return

        self._results = self._parse_results(p3_out, sub_start)
        self._gene_seq = gene_seq
        self._sub_start = sub_start
        self._item_features = getattr(item, "features", [])

        # Populate table
        self.tree.delete(*self.tree.get_children())
        for r in self._results:
            self.tree.insert("", "end", values=(
                r["rank"], r["fwd_seq"], r["rev_seq"],
                r["tm_f"], r["tm_r"], r["gc_f"], r["gc_r"],
                r["size"], r["penalty"]))

        # Auto-select first result
        if self.tree.get_children():
            first = self.tree.get_children()[0]
            self.tree.selection_set(first)
            self._on_result_select()

    # ── Result selection → viewer ─────────────────────────────────────────

    def _on_result_select(self, event=None):
        sel = self.tree.selection()
        if not sel or not self._results:
            return
        idx = self.tree.index(sel[0])
        if idx >= len(self._results):
            return
        r = self._results[idx]
        gene_seq = self._gene_seq

        # Build feature list for viewer
        feats = []
        for feat in self._item_features:
            fs, fe = int(feat.get("start", 0)), int(feat.get("end", 0))
            rel_s = fs - self._sub_start
            rel_e = fe - self._sub_start
            if rel_e > 0 and rel_s < len(gene_seq):
                feats.append({**feat, "start": max(0, rel_s), "end": min(len(gene_seq), rel_e)})

        # Gene target highlight
        feats.append({
            "label": self.feat_c.get() or "Target",
            "start": self._local_start,
            "end": self._local_end,
            "type": "CDS", "strand": 1,
        })

        # Primer positions (relative to gene_seq / sub_start region)
        lpos = r["lpos"]   # [start, len]
        rpos = r["rpos"]   # [3p_end, len]
        fwd_start = lpos[0]
        fwd_end = lpos[0] + lpos[1]
        rev_end = rpos[0] + 1
        rev_start = rpos[0] - rpos[1] + 1

        feats.append({
            "label": f"Fwd#{r['rank']} ({r['tm_f']}°C)",
            "start": fwd_start, "end": fwd_end,
            "type": "Primer", "strand": 1,
        })
        feats.append({
            "label": f"Rev#{r['rank']} ({r['tm_r']}°C)",
            "start": rev_start, "end": rev_end,
            "type": "Primer", "strand": -1,
        })

        UIHelper.render_gene_viewer(self.viewer, gene_seq, feats,
                                    primers_to_show=[], lib_manager=self.lib)

    def _render_viewer(self, feats, primers):
        if self._template_seq:
            UIHelper.render_gene_viewer(self.viewer, self._template_seq[:2000],
                                        feats, primers_to_show=primers,
                                        lib_manager=self.lib)

    # ── Save helpers ──────────────────────────────────────────────────────

    def _save_pair(self, r, prefix):
        fwd_name = f"{prefix}_Fwd"
        rev_name = f"{prefix}_Rev"
        ctr = 2
        base_fwd = fwd_name
        while fwd_name in self.lib.fwd_primers:
            fwd_name = f"{base_fwd} ({ctr})"; ctr += 1
        ctr = 2
        base_rev = rev_name
        while rev_name in self.lib.rev_primers:
            rev_name = f"{base_rev} ({ctr})"; ctr += 1

        fwd_obj = Primer(name=fwd_name, full_sequence=r["fwd_seq"],
                         binding_len=len(r["fwd_seq"]), p_type="Fwd")
        rev_obj = Primer(name=rev_name, full_sequence=r["rev_seq"],
                         binding_len=len(r["rev_seq"]), p_type="Rev")
        self.lib.fwd_primers[fwd_name] = fwd_obj
        self.lib.rev_primers[rev_name] = rev_obj
        return fwd_name, rev_name

    def save_selected(self):
        sel = self.tree.selection()
        if not sel or not self._results:
            messagebox.showwarning("Warning", "저장할 결과를 선택하세요.")
            return
        idx = self.tree.index(sel[0])
        if idx >= len(self._results):
            return
        prefix = self.name_pfx.get().strip() or (self.feat_c.get() or "primer")
        prefix = f"{prefix}_p{idx+1}"
        self.lib.push_undo()
        fwd_n, rev_n = self._save_pair(self._results[idx], prefix)
        self.lib.save()
        self.refresh_cb()
        messagebox.showinfo("Saved", f"저장 완료:\n  Fwd: {fwd_n}\n  Rev: {rev_n}")

    def save_all(self):
        if not self._results:
            messagebox.showwarning("Warning", "먼저 Design Primers를 실행하세요.")
            return
        prefix = self.name_pfx.get().strip() or (self.feat_c.get() or "primer")
        self.lib.push_undo()
        saved = []
        for r in self._results:
            pfx = f"{prefix}_p{r['rank']}"
            fwd_n, rev_n = self._save_pair(r, pfx)
            saved.extend([fwd_n, rev_n])
        self.lib.save()
        self.refresh_cb()
        messagebox.showinfo("Saved", f"{len(saved)}개 primer 저장 완료.")
