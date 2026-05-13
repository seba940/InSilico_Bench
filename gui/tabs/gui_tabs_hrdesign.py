"""HR Primer Auto-Design tab.

Workflow:
  1. Select a genome template with an annotated gene (or enter coordinates manually).
  2. Select a Selection Marker (and optionally an Epitope Tag).
  3. Set homology arm length and marker annealing length.
  4. Click "Design Primers" → get chimeric Fwd / Rev primers ready to save.

Primer structure
----------------
Gene deletion / replacement
  Fwd:  5'─[upstream arm (sense)]─[marker 5' binding]─3'
  Rev:  5'─[RC(downstream arm)]─[RC(marker 3' binding)]─3'

C-terminal tagging (tag inserted before stop codon)
  Fwd:  5'─[before-stop arm (sense)]─[tag 5' binding]─3'
  Rev:  5'─[RC(after-stop arm)]─[RC(tag 3' binding)]─3'
"""
import tkinter as tk
from tkinter import ttk, messagebox
from core.models import Primer, SequenceItem
from core.utils import get_rc
from gui.gui_components import UIHelper, SearchableCombobox, SnapGeneViewer

try:
    from Bio.SeqUtils.MeltingTemp import Tm_NN
except ImportError:
    try:
        from Bio.SeqUtils import Tm_NN
    except ImportError:
        Tm_NN = None


def _calc_tm(seq, binding_len=20):
    if not Tm_NN or not seq:
        return 0.0
    bind = seq[-binding_len:] if len(seq) >= binding_len else seq
    try:
        return round(Tm_NN(bind, Na=50, dnac1=250, dnac2=250), 1)
    except Exception:
        return 0.0


class HRDesignTab:
    """Automatic HR primer designer."""

    def __init__(self, parent, lib, refresh_callback):
        self.parent = parent
        self.lib = lib
        self.refresh_cb = refresh_callback
        self._designed = []   # list of Primer objects from last design run
        self.setup_ui()

    # ── UI setup ─────────────────────────────────────────────────────────

    def setup_ui(self):
        f = ttk.Frame(self.parent, padding=10)
        f.pack(fill="both", expand=True)

        # ── Source section ───────────────────────────────────────────────
        src_f = ttk.LabelFrame(f, text="1. Target Locus", padding=8)
        src_f.pack(fill="x", pady=4)

        r0 = ttk.Frame(src_f); r0.pack(fill="x", pady=2)
        ttk.Label(r0, text="Genome Template:", width=20).pack(side="left")
        self.genome_c = SearchableCombobox(r0)
        self.genome_c.pack(side="left", fill="x", expand=True, padx=4)
        self.genome_c.bind("<<ComboboxSelected>>", lambda e: self._on_genome_change())

        r1 = ttk.Frame(src_f); r1.pack(fill="x", pady=2)
        ttk.Label(r1, text="Annotation Ref:", width=20).pack(side="left")
        self.ann_c = SearchableCombobox(r1)
        self.ann_c.pack(side="left", fill="x", expand=True, padx=4)
        self.ann_c.bind("<<ComboboxSelected>>", lambda e: self._on_ann_change())

        r2 = ttk.Frame(src_f); r2.pack(fill="x", pady=2)
        ttk.Label(r2, text="Gene / Feature:", width=20).pack(side="left")
        self.gene_c = ttk.Combobox(r2, state="readonly")
        self.gene_c.pack(side="left", fill="x", expand=True, padx=4)
        self.gene_c.bind("<<ComboboxSelected>>", lambda e: self._on_gene_change())

        # Manual coordinate override
        r3 = ttk.Frame(src_f); r3.pack(fill="x", pady=2)
        ttk.Label(r3, text="Or manual coords:", width=20).pack(side="left")
        ttk.Label(r3, text="Start:").pack(side="left")
        self.coord_s = ttk.Entry(r3, width=10); self.coord_s.pack(side="left", padx=2)
        ttk.Label(r3, text="End:").pack(side="left")
        self.coord_e = ttk.Entry(r3, width=10); self.coord_e.pack(side="left", padx=2)
        ttk.Label(r3, text="Strand:").pack(side="left")
        self.strand_v = tk.StringVar(value="+")
        ttk.Radiobutton(r3, text="+", variable=self.strand_v, value="+").pack(side="left")
        ttk.Radiobutton(r3, text="-", variable=self.strand_v, value="-").pack(side="left")

        self.locus_info = ttk.Label(src_f, text="", font=("TkDefaultFont", 8, "italic"), foreground="gray")
        self.locus_info.pack(anchor="w")

        # ── Insert section ────────────────────────────────────────────────
        ins_f = ttk.LabelFrame(f, text="2. Insert (Marker / Tag)", padding=8)
        ins_f.pack(fill="x", pady=4)

        ri0 = ttk.Frame(ins_f); ri0.pack(fill="x", pady=2)
        ttk.Label(ri0, text="Design type:", width=20).pack(side="left")
        self.design_type = ttk.Combobox(ri0, state="readonly", width=24,
            values=["Gene Deletion / Replacement",
                    "C-terminal Tagging",
                    "N-terminal Tagging"])
        self.design_type.set("Gene Deletion / Replacement")
        self.design_type.pack(side="left", padx=4)

        ri1 = ttk.Frame(ins_f); ri1.pack(fill="x", pady=2)
        ttk.Label(ri1, text="Selection Marker:", width=20).pack(side="left")
        self.marker_c = SearchableCombobox(ri1)
        self.marker_c.pack(side="left", fill="x", expand=True, padx=4)

        ri2 = ttk.Frame(ins_f); ri2.pack(fill="x", pady=2)
        ttk.Label(ri2, text="Epitope Tag (opt.):", width=20).pack(side="left")
        self.tag_c = SearchableCombobox(ri2)
        self.tag_c.pack(side="left", fill="x", expand=True, padx=4)

        # ── Parameters ──────────────────────────────────────────────────
        par_f = ttk.LabelFrame(f, text="3. Parameters", padding=8)
        par_f.pack(fill="x", pady=4)

        rp = ttk.Frame(par_f); rp.pack(fill="x")
        ttk.Label(rp, text="Homology arm (bp):").pack(side="left")
        self.arm_sp = ttk.Spinbox(rp, from_=20, to=200, width=6)
        self.arm_sp.set(45); self.arm_sp.pack(side="left", padx=4)

        ttk.Label(rp, text="Marker bind (bp):").pack(side="left", padx=(12, 0))
        self.mbind_sp = ttk.Spinbox(rp, from_=15, to=40, width=6)
        self.mbind_sp.set(20); self.mbind_sp.pack(side="left", padx=4)

        ttk.Label(rp, text="Binding_len for Tm:").pack(side="left", padx=(12, 0))
        self.blen_sp = ttk.Spinbox(rp, from_=10, to=40, width=6)
        self.blen_sp.set(20); self.blen_sp.pack(side="left", padx=4)

        # ── Design button ────────────────────────────────────────────────
        ttk.Button(f, text="⚡ Design Primers", command=self.design).pack(pady=6)

        # ── Result section ───────────────────────────────────────────────
        res_f = ttk.LabelFrame(f, text="4. Designed Primers", padding=8)
        res_f.pack(fill="both", expand=True, pady=4)

        cols = ("type", "name", "overhang", "binding", "full", "length", "tm_b", "tm_t")
        heads = ("Type", "Name", "Overhang (5')", "Binding (3')", "Full Sequence", "Len", "Tm(bind)", "Tm(full)")
        fr_t, self.result_tree = UIHelper.create_scrolled_tree(res_f, cols, heads)
        fr_t.pack(fill="x")
        for col, w in zip(cols, [50, 120, 200, 120, 300, 50, 70, 70]):
            self.result_tree.column(col, width=w)

        save_f = ttk.Frame(res_f); save_f.pack(fill="x", pady=4)
        ttk.Label(save_f, text="Name prefix:").pack(side="left")
        self.name_prefix = ttk.Entry(save_f, width=20)
        self.name_prefix.pack(side="left", padx=4)
        ttk.Button(save_f, text="Save Selected to Primer Library",
                   command=self.save_selected).pack(side="left", padx=4)
        ttk.Button(save_f, text="Save All",
                   command=self.save_all).pack(side="left")

        # Viewer to preview the locus with arm annotations
        self.viewer = SnapGeneViewer(res_f, lib_manager=self.lib)
        self.viewer.pack(fill="both", expand=True, pady=4)

    # ── Refresh helpers ───────────────────────────────────────────────────

    def refresh_sources(self):
        def _kind(v):
            return getattr(v, 'kind', None) or ("Genome" if v.topology == "Linear" else "Plasmid")
        gens = sorted(k for k, v in self.lib.templates.items() if _kind(v) == "Genome")
        self.genome_c['values'] = gens
        self.ann_c['values'] = sorted(self.lib.ann_refs.keys())
        self.marker_c['values'] = sorted(self.lib.markers.keys())
        self.tag_c['values'] = [""] + sorted(self.lib.tags.keys())

    def _on_genome_change(self):
        gname = self.genome_c.get()
        if gname in self.lib.templates:
            t = self.lib.templates[gname]
            self.locus_info.config(
                text=f"Genome: {len(t.sequence):,} bp ({t.topology})",
                foreground="black")

    def _on_ann_change(self):
        ann_name = self.ann_c.get()
        if ann_name not in self.lib.ann_refs:
            return
        feats = self.lib.ann_refs[ann_name].features
        gene_names = sorted({f.get("label", "") for f in feats if f.get("label")})
        self.gene_c['values'] = gene_names[:500]

    def _on_gene_change(self):
        ann_name = self.ann_c.get()
        gene_name = self.gene_c.get()
        if not (ann_name in self.lib.ann_refs and gene_name):
            return
        for feat in self.lib.ann_refs[ann_name].features:
            if feat.get("label") == gene_name:
                s, e = feat.get("start", 0), feat.get("end", 0)
                strand = "+" if feat.get("strand", 1) >= 0 else "-"
                self.coord_s.delete(0, tk.END); self.coord_s.insert(0, str(s))
                self.coord_e.delete(0, tk.END); self.coord_e.insert(0, str(e))
                self.strand_v.set(strand)
                self.locus_info.config(
                    text=f"Gene: {gene_name}  {s}–{e} ({strand})  {e-s} bp",
                    foreground="darkgreen")
                break

    # ── Core design logic ─────────────────────────────────────────────────

    def design(self):
        gname = self.genome_c.get()
        marker_name = self.marker_c.get()
        tag_name = self.tag_c.get()
        dtype = self.design_type.get()

        if gname not in self.lib.templates:
            messagebox.showwarning("Warning", "Genome Template을 선택하세요.")
            return
        if not marker_name and dtype != "C-terminal Tagging":
            messagebox.showwarning("Warning", "Selection Marker를 선택하세요.")
            return

        try:
            arm = int(self.arm_sp.get())
            mbind = int(self.mbind_sp.get())
            blen = int(self.blen_sp.get())
            gene_start = int(self.coord_s.get())
            gene_end = int(self.coord_e.get())
        except ValueError:
            messagebox.showerror("Error", "Coordinates와 Parameter가 올바른 숫자인지 확인하세요.")
            return

        genome_seq = self.lib.templates[gname].sequence
        marker_seq = self.lib.markers[marker_name].sequence if marker_name in self.lib.markers else ""
        tag_seq = self.lib.tags[tag_name].sequence if tag_name in self.lib.tags else ""

        on_minus = (self.strand_v.get() == "-")

        # Clamp to genome boundaries
        gene_start = max(0, gene_start)
        gene_end = min(len(genome_seq), gene_end)

        if dtype == "Gene Deletion / Replacement":
            insert_seq = marker_seq
            # Fwd: upstream arm + marker 5' bind
            upstream = genome_seq[max(0, gene_start - arm): gene_start]
            fwd_full = upstream + insert_seq[:mbind]
            # Rev: RC(downstream arm) + RC(marker 3' bind)
            downstream = genome_seq[gene_end: gene_end + arm]
            rev_full = get_rc(downstream) + get_rc(insert_seq[-mbind:])
            if on_minus:
                fwd_full, rev_full = get_rc(rev_full), get_rc(fwd_full)

        elif dtype == "C-terminal Tagging":
            insert_seq = tag_seq if tag_seq else marker_seq
            if not insert_seq:
                messagebox.showwarning("Warning", "Epitope Tag 또는 Marker를 선택하세요.")
                return
            # Insert tag just before stop codon
            # Fwd arm: sequence just upstream of stop (within the gene, near end)
            before_stop = genome_seq[max(0, gene_end - 3 - arm): gene_end - 3]
            fwd_full = before_stop + insert_seq[:mbind]
            # Rev arm: RC of sequence just after stop codon
            after_stop = genome_seq[gene_end: gene_end + arm]
            rev_full = get_rc(after_stop) + get_rc(insert_seq[-mbind:])
            if on_minus:
                fwd_full, rev_full = get_rc(rev_full), get_rc(fwd_full)

        elif dtype == "N-terminal Tagging":
            insert_seq = tag_seq if tag_seq else marker_seq
            if not insert_seq:
                messagebox.showwarning("Warning", "Epitope Tag 또는 Marker를 선택하세요.")
                return
            # Insert tag just after start codon (ATG)
            after_atg = genome_seq[gene_start + 3: gene_start + 3 + arm]
            fwd_full = after_atg + insert_seq[:mbind]
            before_atg = genome_seq[max(0, gene_start - arm): gene_start]
            rev_full = get_rc(before_atg) + get_rc(insert_seq[-mbind:])
            if on_minus:
                fwd_full, rev_full = get_rc(rev_full), get_rc(fwd_full)
        else:
            messagebox.showerror("Error", f"Unknown design type: {dtype}")
            return

        prefix = self.name_prefix.get().strip() or (self.gene_c.get() or "HR")
        fwd_name = f"{prefix}_Fwd"
        rev_name = f"{prefix}_Rev"

        fwd_obj = Primer(
            name=fwd_name,
            full_sequence=fwd_full.upper(),
            binding_len=min(blen, len(fwd_full)),
            p_type="Fwd")
        rev_obj = Primer(
            name=rev_name,
            full_sequence=rev_full.upper(),
            binding_len=min(blen, len(rev_full)),
            p_type="Rev")

        self._designed = [fwd_obj, rev_obj]

        # Populate result tree
        self.result_tree.delete(*self.result_tree.get_children())
        for p in [fwd_obj, rev_obj]:
            self.result_tree.insert("", "end", values=(
                p.p_type, p.name,
                p.overhang, p.binding,
                p.full_sequence, len(p.full_sequence),
                p.tm, p.tm_total))

        # Preview locus in viewer
        ext_s = max(0, gene_start - arm - 20)
        ext_e = min(len(genome_seq), gene_end + arm + 20)
        ext_seq = genome_seq[ext_s:ext_e]
        feats = [
            {"label": self.gene_c.get() or "Target",
             "start": gene_start - ext_s,
             "end": gene_end - ext_s,
             "type": "CDS", "strand": 1 if not on_minus else -1},
            {"label": f"Fwd arm ({arm}bp)",
             "start": max(0, gene_start - arm) - ext_s,
             "end": gene_start - ext_s,
             "type": "Homology Arm", "strand": 1},
            {"label": f"Rev arm ({arm}bp)",
             "start": gene_end - ext_s,
             "end": min(len(ext_seq), gene_end + arm - ext_s),
             "type": "Homology Arm", "strand": -1},
        ]
        UIHelper.render_gene_viewer(self.viewer, ext_seq, feats,
                                    primers_to_show=[fwd_obj, rev_obj],
                                    lib_manager=self.lib)

    # ── Save helpers ─────────────────────────────────────────────────────

    def _save_primer(self, p):
        name = p.name; ctr = 2
        bucket = self.lib.fwd_primers if p.p_type == "Fwd" else self.lib.rev_primers
        while name in bucket:
            name = f"{p.name} ({ctr})"; ctr += 1
        p.name = name
        bucket[name] = p

    def save_selected(self):
        if not self._designed:
            messagebox.showwarning("Warning", "먼저 Design Primers를 실행하세요.")
            return
        sel = self.result_tree.selection()
        indices = [self.result_tree.index(s) for s in sel]
        primers = [self._designed[i] for i in indices if i < len(self._designed)]
        if not primers:
            messagebox.showwarning("Warning", "저장할 프라이머를 선택하세요.")
            return
        self.lib.push_undo()
        for p in primers:
            self._save_primer(p)
        self.lib.save()
        self.refresh_cb()
        messagebox.showinfo("Saved", f"{len(primers)}개 primer가 라이브러리에 저장되었습니다.")

    def save_all(self):
        if not self._designed:
            messagebox.showwarning("Warning", "먼저 Design Primers를 실행하세요.")
            return
        self.lib.push_undo()
        for p in self._designed:
            self._save_primer(p)
        self.lib.save()
        self.refresh_cb()
        messagebox.showinfo("Saved", f"{len(self._designed)}개 primer가 라이브러리에 저장되었습니다.")
