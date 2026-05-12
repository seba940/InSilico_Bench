"""Gene Lookup tab - SGD/NCBI fetch"""
import webbrowser
import tkinter as tk
from tkinter import ttk, messagebox
from core.models import SequenceItem
from gui.gui_components import UIHelper, SearchableCombobox, SnapGeneViewer
from core.utils import get_rc

SGD_LOCUS = "https://www.yeastgenome.org/locus/"

class GeneLookupTab:
    def __init__(self, parent, lib, refresh_callback):
        self.parent = parent
        self.lib = lib
        self.refresh_cb = refresh_callback
        self.matches = []
        self.current_extract = None
        self.setup_ui()

    def setup_ui(self):
        f = ttk.Frame(self.parent, padding=10)
        f.pack(fill="both", expand=True)
        sf = ttk.LabelFrame(f, text="Source", padding=8)
        sf.pack(fill="x")
        r0 = ttk.Frame(sf); r0.pack(fill="x", pady=2)
        ttk.Label(r0, text="Genome template:", width=18).pack(side="left")
        self.gen_c = SearchableCombobox(r0)
        self.gen_c.pack(side="left", fill="x", expand=True, padx=4)
        self.gen_c.bind("<<ComboboxSelected>>", lambda e: self._update_status())
        r1 = ttk.Frame(sf); r1.pack(fill="x", pady=2)
        ttk.Label(r1, text="Annotation reference:", width=18).pack(side="left")
        self.ann_c = SearchableCombobox(r1)
        self.ann_c.pack(side="left", fill="x", expand=True, padx=4)
        self.ann_c.bind("<<ComboboxSelected>>", lambda e: self.refresh_results())
        self.status_l = ttk.Label(sf, text="", font=("TkDefaultFont", 8, "italic"), foreground="gray")
        self.status_l.pack(anchor="w", pady=2)
        srch = ttk.LabelFrame(f, text="Search Gene", padding=8)
        srch.pack(fill="x", pady=4)
        rs = ttk.Frame(srch); rs.pack(fill="x")
        ttk.Label(rs, text="Name / Alias contains:", width=22).pack(side="left")
        self.q = ttk.Entry(rs)
        self.q.pack(side="left", fill="x", expand=True, padx=4)
        self.q.bind("<KeyRelease>", lambda e: self.refresh_results())
        ttk.Label(rs, text="Type:").pack(side="left", padx=(8, 0))
        self.type_filter = ttk.Combobox(rs, state="readonly", width=10,
            values=["(all)","gene","CDS","tRNA","rRNA","ncRNA","snoRNA","telomere","centromere","ARS"])
        self.type_filter.set("gene")
        self.type_filter.pack(side="left", padx=4)
        self.type_filter.bind("<<ComboboxSelected>>", lambda e: self.refresh_results())
        ft, self.tree = UIHelper.create_scrolled_tree(f,
            ("name","sysname","type","chr","s","e","strand","len"),
            ("Name","Sysname","Type","Chr","Start","End","Strand","Length"))
        ft.pack(fill="x", pady=4)
        for col, w in [("name",100),("sysname",100),("type",70),("chr",70),
                       ("s",80),("e",80),("strand",60),("len",80)]:
            self.tree.column(col, width=w)
        self.tree.bind("<<TreeviewSelect>>", lambda e: self.on_select())
        ff = ttk.LabelFrame(f, text="Extract with Flanking", padding=8)
        ff.pack(fill="x", pady=4)
        rf = ttk.Frame(ff); rf.pack(fill="x")
        ttk.Label(rf, text="Upstream:").pack(side="left")
        self.up_sp = ttk.Spinbox(rf, from_=0, to=20000, width=8, command=self.on_select)
        self.up_sp.set(500); self.up_sp.pack(side="left", padx=4)
        self.up_sp.bind("<Return>", lambda e: self.on_select())
        ttk.Label(rf, text="bp  Downstream:").pack(side="left")
        self.dn_sp = ttk.Spinbox(rf, from_=0, to=20000, width=8, command=self.on_select)
        self.dn_sp.set(500); self.dn_sp.pack(side="left", padx=4)
        self.dn_sp.bind("<Return>", lambda e: self.on_select())
        ttk.Label(rf, text="bp").pack(side="left")
        self.match_l = ttk.Label(rf, text="", font=("TkDefaultFont",9,"bold"), foreground="darkgreen")
        self.match_l.pack(side="right", padx=10)
        self.det = SnapGeneViewer(ff, row_width=60)
        self.det.pack(fill="both", expand=True, pady=4)
        sgd_f = ttk.LabelFrame(f, text="SGD Database", padding=8)
        sgd_f.pack(fill="x", pady=4)
        rs2 = ttk.Frame(sgd_f); rs2.pack(fill="x")
        ttk.Button(rs2, text="Open in SGD", command=self.open_in_sgd).pack(side="left", padx=4)
        ttk.Button(rs2, text="Open SGD Sequence Page", command=self.open_sgd_seq).pack(side="left", padx=4)
        self.sgd_status = ttk.Label(rs2, text="", font=("TkDefaultFont",9), foreground="gray")
        self.sgd_status.pack(side="left", padx=8)
        sv = ttk.Frame(f); sv.pack(fill="x", pady=4)
        ttk.Label(sv, text="Save name:").pack(side="left")
        self.save_name = ttk.Entry(sv, width=40)
        self.save_name.pack(side="left", padx=4)
        ttk.Button(sv, text="Save as Template (Kind=Genome)",
                   command=self.save_as_template).pack(side="left", padx=4)
        # legend is now embedded inside SnapGeneViewer toolbar

    def refresh_sources(self):
        def _kind(v):
            return v.kind or ("Genome" if v.topology == "Linear" else "Plasmid")
        gens = [k for k,v in self.lib.templates.items() if _kind(v) == "Genome"]
        self.gen_c["values"] = sorted(gens)
        self.ann_c["values"] = sorted(list(self.lib.ann_refs.keys()))

    def _update_status(self):
        gname = self.gen_c.get()
        if gname and gname in self.lib.templates:
            t = self.lib.templates[gname]
            self.status_l.config(
                text="Genome: {:,} bp ({})".format(len(t.sequence), t.topology),
                foreground="black")
        else:
            self.status_l.config(text="", foreground="gray")
        self.refresh_results()

    def refresh_results(self):
        self.tree.delete(*self.tree.get_children())
        ann_name = self.ann_c.get()
        if not ann_name or ann_name not in self.lib.ann_refs:
            return
        feats = self.lib.ann_refs[ann_name].features
        q = self.q.get().strip().upper()
        tfilter = self.type_filter.get()
        self.matches = []
        for feat in feats:
            if tfilter != "(all)" and feat.get("type") != tfilter:
                continue
            label = (feat.get("label") or "").upper()
            aliases_up = [a.upper() for a in feat.get("aliases", [])]
            if q and q not in " ".join([label] + aliases_up):
                continue
            self.matches.append(feat)
        cap = 500
        for i, feat in enumerate(self.matches[:cap]):
            strand = "+" if feat.get("strand", 1) == 1 else "-"
            length = feat.get("end", 0) - feat.get("start", 0)
            attrs = feat.get("attrs", {}) or {}
            sysname = attrs.get("Name") or attrs.get("ID") or ""
            self.tree.insert("", "end", iid=str(i), values=(
                feat.get("label",""), sysname, feat.get("type",""),
                feat.get("chr",""), feat.get("start",""), feat.get("end",""),
                strand, length))
        count = len(self.matches)
        if count > cap:
            txt = "{} hits (showing first {})".format(count, cap)
        else:
            txt = "{} hit(s)".format(count)
        self.match_l.config(text=txt, foreground="darkgreen")

    def on_select(self):
        sel = self.tree.selection()
        if not sel:
            return
        try:
            i = int(sel[0])
        except ValueError:
            return
        if i >= len(self.matches):
            return
        feat = self.matches[i]
        gname = self.gen_c.get()
        if not gname or gname not in self.lib.templates:
            self.det.show_message("Genome template을 먼저 선택하세요.")
            return
        gt = self.lib.templates[gname]
        gseq_up = gt.sequence.upper()
        gene_seq = (feat.get("sequence") or "").upper()
        try:
            up = max(0, int(self.up_sp.get()))
            dn = max(0, int(self.dn_sp.get()))
        except ValueError:
            up, dn = 500, 500
        gene_start, gene_end, strand_in_genome = -1, -1, +1
        if gene_seq:
            idx = gseq_up.find(gene_seq)
            if idx != -1:
                gene_start, gene_end, strand_in_genome = idx, idx + len(gene_seq), +1
            else:
                rc = get_rc(gene_seq).upper()
                idx_rc = gseq_up.find(rc)
                if idx_rc != -1:
                    gene_start, gene_end, strand_in_genome = idx_rc, idx_rc + len(rc), -1
        if gene_start == -1:
            f_start = feat.get("start")
            f_end = feat.get("end")
            if f_start is not None and f_end is not None and f_end > f_start:
                gene_start = max(0, f_start)
                gene_end = min(len(gt.sequence), f_end)
                strand_in_genome = feat.get("strand", 1)
        if gene_start == -1:
            _msg = ("genome '{gn}' 에서 '{lb}' 을 찾지 못했습니다.\n"
                   "GFF 좌표: {ch}:{s}..{e}\n"
                   "다른 strain 이거나 해당 chromosome 을 포함하지 않음.").format(
                       gn=gname, lb=feat.get("label"),
                       ch=feat.get("chr"), s=feat.get("start"), e=feat.get("end"))
            self.det.show_message(_msg)
            self.match_l.config(text="not found in genome", foreground="red")
            self.current_extract = None
            return
        ext_start = max(0, gene_start - up)
        ext_end = min(len(gt.sequence), gene_end + dn)
        ext_seq = gt.sequence[ext_start:ext_end]
        gls = gene_start - ext_start
        gle = gene_end - ext_start
        ft2 = feat.get("type","CDS")
        if ft2 not in ("gene","CDS"):
            ft2 = "CDS"
        feats_local = [{"label": feat.get("label","gene"),
                        "start": gls, "end": gle, "type": ft2,
                        "strand": strand_in_genome}]
        if up > 0:
            feats_local.append({"label": "upstream_{}bp".format(up),
                                 "start": 0, "end": gls, "type": "Misc", "strand": 1})
        if dn > 0:
            feats_local.append({"label": "downstream_{}bp".format(dn),
                                 "start": gle, "end": len(ext_seq), "type": "Misc", "strand": 1})
        UIHelper.render_annotations(self.det, ext_seq, feats_local, lib_manager=self.lib)
        ss = "+" if strand_in_genome == +1 else "-"
        self.match_l.config(
            text="matched {}..{} ({}), extract {} bp".format(gene_start, gene_end, ss, len(ext_seq)),
            foreground="darkgreen")
        self.current_extract = {
            "name_default": "{}_{}up_{}dn".format(feat.get("label","gene"), up, dn),
            "sequence": ext_seq, "features": feats_local,
            "src_genome": gname, "gene_label": feat.get("label","gene")}
        self.save_name.delete(0, tk.END)
        self.save_name.insert(0, self.current_extract["name_default"])

    # -------- SGD 브라우저 열기 ----------------------------------------
    def _get_sgd_identifier(self):
        """선택된 유전자의 SGD 식별자 반환. 없으면 None."""
        sel = self.tree.selection()
        if not sel:
            messagebox.showwarning("Warning", "유전자를 먼저 선택하세요.")
            return None
        try:
            i = int(sel[0])
        except ValueError:
            return None
        if i >= len(self.matches):
            return None
        feat = self.matches[i]
        attrs = feat.get("attrs", {}) or {}
        # standard name(ACT1) 우선, 없으면 systematic(YFL039C)
        ident = (attrs.get("gene") or feat.get("label","") or
                 attrs.get("Name") or attrs.get("ID") or "").strip()
        return ident or None

    def open_in_sgd(self):
        """선택된 유전자의 SGD locus 페이지를 기본 브라우저로 연다."""
        ident = self._get_sgd_identifier()
        if not ident:
            return
        url = SGD_LOCUS + ident
        webbrowser.open(url)
        self.sgd_status.config(text="Opened: " + url, foreground="darkgreen")

    def open_sgd_seq(self):
        """SGD locus 페이지의 Sequence 탭을 바로 연다."""
        ident = self._get_sgd_identifier()
        if not ident:
            return
        url = SGD_LOCUS + ident + "#sequence"
        webbrowser.open(url)
        self.sgd_status.config(text="Opened: " + url, foreground="darkgreen")

    def save_as_template(self):
        if not self.current_extract:
            messagebox.showwarning("Warning", "추출할 유전자를 먼저 선택하세요.")
            return
        name = self.save_name.get().strip() or self.current_extract["name_default"]
        original = name; ctr = 2
        while name in self.lib.templates:
            name = "{} ({})".format(original, ctr); ctr += 1
        src = self.lib.templates.get(self.current_extract["src_genome"])
        species = src.species if src else ""
        item = SequenceItem(
            name=name, sequence=self.current_extract["sequence"],
            category="template", topology="Linear", kind="Genome",
            species=species, features=self.current_extract["features"],
            template_name=self.current_extract["src_genome"])
        self.lib.templates[name] = item
        self.lib.save()
        self.refresh_cb()
        msg = "'" + name + "' Template으로 저장 완료."
        messagebox.showinfo("Saved", msg)
