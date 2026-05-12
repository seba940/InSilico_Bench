import tkinter as tk
from tkinter import ttk, messagebox, scrolledtext, filedialog
from core.models import Primer, SequenceItem, Marker, EpitopeTag, GlobalFeature
from gui.gui_components import UIHelper, FeatureEditor, SearchableCombobox, FindDialog, SnapGeneViewer
import pandas as pd
from Bio import SeqIO
import os

class TemplateTab:
    def __init__(self, parent, lib, refresh_callback):
        self.parent, self.lib, self.refresh_cb = parent, lib, refresh_callback
        self.setup_ui()

    def setup_ui(self):
        f = ttk.Frame(self.parent, padding=10); f.pack(fill="both", expand=True)
        pane = ttk.PanedWindow(f, orient="horizontal")
        pane.pack(fill="both", expand=True)
        l_f = ttk.Frame(pane)
        pane.add(l_f, weight=1)
        r_f = ttk.Frame(pane)
        pane.add(r_f, weight=2)
        in_f = ttk.LabelFrame(l_f, text="Register Template", padding=10); in_f.pack(fill="x")
        ttk.Label(in_f, text="Name:").pack(anchor="w"); self.tn = ttk.Entry(in_f); self.tn.pack(fill="x")
        self.ts = scrolledtext.ScrolledText(in_f, height=5, exportselection=False); self.ts.pack(fill="x")
        # Topology
        r1 = ttk.Frame(in_f); r1.pack(fill="x", pady=2)
        ttk.Label(r1, text="Topology:", width=10).pack(side="left")
        self.topo_v = tk.StringVar(value="Linear")
        ttk.Radiobutton(r1, text="Linear", variable=self.topo_v, value="Linear",
                        command=self._sync_kind_from_topo).pack(side="left")
        ttk.Radiobutton(r1, text="Circular", variable=self.topo_v, value="Circular",
                        command=self._sync_kind_from_topo).pack(side="left")
        # Kind (Genome/Plasmid) — HR 탭에서 target 분류에 사용
        r1b = ttk.Frame(in_f); r1b.pack(fill="x", pady=2)
        ttk.Label(r1b, text="Kind:", width=10).pack(side="left")
        self.kind_v = tk.StringVar(value="Genome")
        ttk.Radiobutton(r1b, text="Genome (HR target)", variable=self.kind_v, value="Genome").pack(side="left")
        ttk.Radiobutton(r1b, text="Plasmid (PCR/Digest source)", variable=self.kind_v, value="Plasmid").pack(side="left")
        # Species
        r2 = ttk.Frame(in_f); r2.pack(fill="x", pady=2)
        ttk.Label(r2, text="Species:", width=10).pack(side="left")
        self.sp_c = SearchableCombobox(r2); self.sp_c.pack(side="left", fill="x", expand=True)
        btn_row = ttk.Frame(in_f); btn_row.pack(fill="x", pady=(4, 0))
        ttk.Button(btn_row, text="Import File (.fasta/.fa/.fsa/.dna)",
                   command=self.import_from_file).pack(side="left")
        ttk.Button(btn_row, text="Save Template", command=self.add_t).pack(side="right")
        self.det = SnapGeneViewer(r_f); self.det.pack(fill="both", expand=True)

        self.feat_ed = FeatureEditor(l_f, "Local Features  (Ctrl+F: 서열 검색 → 위치 자동입력)", text_widget=self.ts); self.feat_ed.pack(fill="x")
        self.ts.bind("<Control-f>", lambda e: FindDialog(self.parent, self.ts, self.feat_ed) or "break")
        self.ts.bind("<Control-F>", lambda e: FindDialog(self.parent, self.ts, self.feat_ed) or "break")
        self.search = UIHelper.add_search_bar(l_f, self.filter, lambda: self.remove_dups())
        fr_t, self.tree = UIHelper.create_scrolled_tree(
            l_f, ("N", "L", "T", "K", "S"),
            ("Name", "Length", "Topology", "Kind", "Species"),
            lib=self.lib, table_id="template_main")
        fr_t.pack(fill="both", expand=True); self.tree.bind("<<TreeviewSelect>>", self.on_sel)
        bf = ttk.Frame(l_f); bf.pack(fill="x")
        ttk.Button(bf, text="Edit Selected", command=self.edit_load).pack(side="left"); ttk.Button(bf, text="Delete", command=self.del_t).pack(side="right")

    def _sync_kind_from_topo(self):
        """Topology 변경 시 Kind 기본값을 자연스럽게 맞춰줌 (사용자가 다시 바꿀 수 있음)."""
        self.kind_v.set("Genome" if self.topo_v.get() == "Linear" else "Plasmid")

    def edit_load(self):
        sel = self.tree.selection()
        if sel:
            i = self.lib.templates[self.tree.item(sel[0])['values'][0]]; self.tn.delete(0, tk.END); self.tn.insert(0, i.name)
            self.ts.delete("1.0", tk.END); self.ts.insert(tk.END, i.sequence); self.topo_v.set(i.topology); self.sp_c.set(i.species); self.feat_ed.set_features(i.features)
            self.kind_v.set(i.kind or ("Genome" if i.topology == "Linear" else "Plasmid"))

    def import_from_file(self):
        """FASTA (.fasta/.fa/.fsa) 또는 SnapGene (.dna) 파일을 열어 서열을 가져온다.
        단일 레코드면 Name/Sequence 필드에 채워 주고,
        다중 레코드(FASTA)면 모두 라이브러리에 바로 저장한다."""
        path = filedialog.askopenfilename(
            title="Open sequence file",
            filetypes=[
                ("Sequence files", "*.fasta *.fa *.fsa *.dna"),
                ("FASTA", "*.fasta *.fa *.fsa"),
                ("SnapGene", "*.dna"),
                ("All files", "*.*"),
            ])
        if not path:
            return

        ext = os.path.splitext(path)[1].lower()
        fmt = "snapgene" if ext == ".dna" else "fasta"

        try:
            records = list(SeqIO.parse(path, fmt))
        except Exception as e:
            messagebox.showerror("Import Error", "파일을 읽지 못했습니다:\n" + str(e))
            return

        if not records:
            messagebox.showwarning("Import", "서열 레코드를 찾지 못했습니다.")
            return

        def _seq_features(record):
            """BioPython SeqFeature 목록 → 앱 feature dict 목록으로 변환."""
            result = []
            for f in getattr(record, "features", []):
                label = (f.qualifiers.get("label") or
                         f.qualifiers.get("gene") or
                         f.qualifiers.get("product") or [f.type])
                label = label[0] if isinstance(label, list) else label
                start = int(f.location.start)
                end = int(f.location.end)
                strand = f.location.strand if f.location.strand is not None else 1
                ftype = f.type if f.type in (
                    "gene", "CDS", "promoter", "terminator",
                    "rep_origin", "misc_feature", "RBS") else "Misc"
                if end > start:
                    result.append({"label": label, "start": start,
                                   "end": end, "type": ftype, "strand": strand})
            return result

        def _topology(record):
            ann = getattr(record, "annotations", {})
            mol = ann.get("molecule_type", "")
            topology = ann.get("topology", "")
            if topology.lower() == "circular" or "circular" in mol.lower():
                return "Circular"
            return "Linear"

        if len(records) == 1:
            # 단일 레코드: 입력 필드에 채워주기
            rec = records[0]
            seq = str(rec.seq).upper()
            name = rec.name or rec.id or os.path.splitext(os.path.basename(path))[0]
            topo = _topology(rec)
            feats = _seq_features(rec)

            self.tn.delete(0, tk.END)
            self.tn.insert(0, name)
            self.ts.delete("1.0", tk.END)
            self.ts.insert(tk.END, seq)
            self.topo_v.set(topo)
            self.kind_v.set("Genome" if topo == "Linear" else "Plasmid")
            self.feat_ed.set_features(feats)
            messagebox.showinfo(
                "Import OK",
                "'{}'  ({:,} bp, {})\n\n"
                "Name/Sequence 필드에 불러왔습니다.\n"
                "'Save Template' 버튼으로 저장하세요.".format(name, len(seq), topo))
        else:
            # 다중 레코드: 전부 라이브러리에 저장
            saved, skipped = 0, 0
            for rec in records:
                seq = str(rec.seq).upper()
                if not seq:
                    skipped += 1
                    continue
                name = rec.name or rec.id or "seq_{}".format(saved + 1)
                # 중복 이름 처리
                orig = name; ctr = 2
                while name in self.lib.templates:
                    name = "{} ({})".format(orig, ctr); ctr += 1
                topo = _topology(rec)
                feats = _seq_features(rec)
                self.lib.templates[name] = SequenceItem(
                    name, seq,
                    topology=topo,
                    kind="Genome" if topo == "Linear" else "Plasmid",
                    features=feats,
                )
                saved += 1
            self.lib.save()
            self.refresh_cb()
            messagebox.showinfo(
                "Import OK",
                "{} 개 레코드를 Template 라이브러리에 저장했습니다.{}".format(
                    saved,
                    "\n({} 개 빈 서열 건너뜀)".format(skipped) if skipped else ""))

    def add_t(self):
        n, s = self.tn.get().strip(), self.ts.get("1.0", tk.END).strip()
        if n and s:
            self.lib.templates[n] = SequenceItem(
                n, s, topology=self.topo_v.get(), species=self.sp_c.get(),
                features=self.feat_ed.get_features(), kind=self.kind_v.get(),
            )
            self.lib.save(); self.refresh_cb()

    def del_t(self):
        for s in self.tree.selection(): del self.lib.templates[self.tree.item(s)['values'][0]]
        self.lib.save(); self.refresh_cb()

    def on_sel(self, e):
        sel = self.tree.selection()
        if sel: 
            i = self.lib.templates[self.tree.item(sel[0])['values'][0]]
            UIHelper.render_annotations(self.det, i.sequence, i.features, lib_manager=self.lib)
            self.feat_ed.set_features(i.features)

    def filter(self, q):
        self.tree.delete(*self.tree.get_children())
        for k, v in self.lib.templates.items():
            if not q or q.upper() in k.upper():
                kind = v.kind or ("Genome" if v.topology == "Linear" else "Plasmid")
                self.tree.insert("", "end", values=(v.name, len(v.sequence), v.topology, kind, v.species))

    def remove_dups(self):
        seen, dups = set(), []
        for k, v in self.lib.templates.items():
            if v.sequence in seen: dups.append(k)
            else: seen.add(v.sequence)
        for k in dups: del self.lib.templates[k]
        if dups: self.lib.save(); self.refresh_cb()

class MarkerTab:
    def __init__(self, parent, lib, refresh_callback):
        self.parent, self.lib, self.refresh_cb = parent, lib, refresh_callback
        self.setup_ui()

    def setup_ui(self):
        f = ttk.Frame(self.parent, padding=10); f.pack(fill="both", expand=True)
        pane = ttk.PanedWindow(f, orient="horizontal")
        pane.pack(fill="both", expand=True)
        l = ttk.Frame(pane)
        pane.add(l, weight=1)
        r_f = ttk.Frame(pane)
        pane.add(r_f, weight=2)
        in_f = ttk.LabelFrame(l, text="Marker Registration", padding=10); in_f.pack(fill="x")
        ttk.Label(in_f, text="Name:").pack(anchor="w"); self.mn = ttk.Entry(in_f); self.mn.pack(fill="x")
        ttk.Label(in_f, text="Sequence:").pack(anchor="w"); self.ms = scrolledtext.ScrolledText(in_f, height=4, exportselection=False); self.ms.pack(fill="x")
        self.det = SnapGeneViewer(r_f); self.det.pack(fill="both", expand=True)
        self.feat_ed = FeatureEditor(l, "Internal Marker Features", text_widget=self.ms); self.feat_ed.pack(fill="x")
        ttk.Button(l, text="Save Marker to Global Library", command=self.add_m).pack()
        self.search = UIHelper.add_search_bar(l, self.filter, self.remove_dups)
        fr_t, self.tree = UIHelper.create_scrolled_tree(l, ("N", "L"), ("Name", "Length"), lib=self.lib, table_id="marker_or_tag")
        fr_t.pack(fill="both", expand=True); self.tree.bind("<<TreeviewSelect>>", self.on_sel)
        ttk.Button(l, text="Edit Selected", command=self.edit_load).pack(side="left"); ttk.Button(l, text="Delete", command=self.del_m).pack(side="right")

    def edit_load(self):
        sel = self.tree.selection()
        if sel:
            m = self.lib.markers[self.tree.item(sel[0])['values'][0]]
            self.mn.delete(0, tk.END); self.mn.insert(0, m.name); self.ms.delete("1.0", tk.END); self.ms.insert(tk.END, m.sequence); self.feat_ed.set_features(m.features)

    def add_m(self):
        n, s = self.mn.get().strip(), self.ms.get("1.0", tk.END).strip()
        if n and s: 
            self.lib.markers[n] = Marker(n, s, features=self.feat_ed.get_features())
            # 중앙 라이브러리 자동 등록
            self.lib.global_features[f"Marker_{n}"] = GlobalFeature(n, s, "Marker")
            self.lib.save(); self.refresh_cb()

    def del_m(self):
        for s in self.tree.selection(): 
            name = self.tree.item(s)['values'][0]
            del self.lib.markers[name]
            if f"Marker_{name}" in self.lib.global_features: del self.lib.global_features[f"Marker_{name}"]
        self.lib.save(); self.refresh_cb()

    def on_sel(self, e):
        sel = self.tree.selection()
        if sel: 
            m = self.lib.markers[self.tree.item(sel[0])['values'][0]]
            UIHelper.render_annotations(self.det, m.sequence, m.features, lib_manager=None)
            self.feat_ed.set_features(m.features)

    def filter(self, q):
        self.tree.delete(*self.tree.get_children())
        for k, v in self.lib.markers.items():
            if not q or q.upper() in k.upper(): self.tree.insert("", "end", values=(v.name, len(v.sequence)))

    def remove_dups(self):
        seen, dups = set(), []
        for k, v in self.lib.markers.items():
            if v.sequence in seen: dups.append(k)
            else: seen.add(v.sequence)
        for k in dups: del self.lib.markers[k]
        if dups: self.lib.save(); self.refresh_cb()

class TagTab:
    def __init__(self, parent, lib, refresh_callback):
        self.parent, self.lib, self.refresh_cb = parent, lib, refresh_callback
        self.setup_ui()

    def setup_ui(self):
        f = ttk.Frame(self.parent, padding=10); f.pack(fill="both", expand=True)
        pane = ttk.PanedWindow(f, orient="horizontal")
        pane.pack(fill="both", expand=True)
        l = ttk.Frame(pane)
        pane.add(l, weight=1)
        r_f = ttk.Frame(pane)
        pane.add(r_f, weight=2)
        in_f = ttk.LabelFrame(l, text="Tag Registration", padding=10); in_f.pack(fill="x")
        ttk.Label(in_f, text="Name:").pack(anchor="w"); self.tn = ttk.Entry(in_f); self.tn.pack(fill="x")
        ttk.Label(in_f, text="Sequence:").pack(anchor="w"); self.ts = scrolledtext.ScrolledText(in_f, height=4, exportselection=False); self.ts.pack(fill="x")
        self.det = SnapGeneViewer(r_f); self.det.pack(fill="both", expand=True)
        self.feat_ed = FeatureEditor(l, "Internal Tag Features", text_widget=self.ts); self.feat_ed.pack(fill="x")
        ttk.Button(l, text="Save Tag to Global Library", command=self.add_t).pack()
        self.search = UIHelper.add_search_bar(l, self.filter, self.remove_dups)
        fr_t, self.tree = UIHelper.create_scrolled_tree(l, ("N", "L"), ("Name", "Length"), lib=self.lib, table_id="marker_or_tag")
        fr_t.pack(fill="both", expand=True); self.tree.bind("<<TreeviewSelect>>", self.on_sel)
        ttk.Button(l, text="Edit Selected", command=self.edit_load).pack(side="left"); ttk.Button(l, text="Delete", command=self.del_t).pack(side="right")

    def edit_load(self):
        sel = self.tree.selection()
        if sel:
            t = self.lib.tags[self.tree.item(sel[0])['values'][0]]
            self.tn.delete(0, tk.END); self.tn.insert(0, t.name); self.ts.delete("1.0", tk.END); self.ts.insert(tk.END, t.sequence); self.feat_ed.set_features(t.features)

    def add_t(self):
        n, s = self.tn.get().strip(), self.ts.get("1.0", tk.END).strip()
        if n and s: 
            self.lib.tags[n] = EpitopeTag(n, s, features=self.feat_ed.get_features())
            self.lib.global_features[f"Tag_{n}"] = GlobalFeature(n, s, "Tag")
            self.lib.save(); self.refresh_cb()

    def del_t(self):
        for s in self.tree.selection(): 
            name = self.tree.item(s)['values'][0]
            del self.lib.tags[name]
            if f"Tag_{name}" in self.lib.global_features: del self.lib.global_features[f"Tag_{name}"]
        self.lib.save(); self.refresh_cb()

    def on_sel(self, e):
        sel = self.tree.selection()
        if sel: 
            t = self.lib.tags[self.tree.item(sel[0])['values'][0]]
            UIHelper.render_annotations(self.det, t.sequence, t.features, lib_manager=None)
            self.feat_ed.set_features(t.features)

    def filter(self, q):
        self.tree.delete(*self.tree.get_children())
        for k, v in self.lib.tags.items():
            if not q or q.upper() in k.upper(): self.tree.insert("", "end", values=(v.name, len(v.sequence)))

    def remove_dups(self):
        seen, dups = set(), []
        for k, v in self.lib.tags.items():
            if v.sequence in seen: dups.append(k)
            else: seen.add(v.sequence)
        for k in dups: del self.lib.tags[k]
        if dups: self.lib.save(); self.refresh_cb()

from core.utils import get_primer_analysis

class PrimerTab:
    def __init__(self, parent, lib, refresh_callback):
        self.parent, self.lib, self.refresh_cb, self.edit_widgets = parent, lib, refresh_callback, []
        self.validated = False  # Validation status
        self.last_validated_seq = ""
        self.setup_ui()

    def setup_ui(self):
        f = ttk.Frame(self.parent, padding=10); f.pack(fill="both", expand=True)
        # Tm 계산 안내
        ff = ttk.LabelFrame(f, text="Tm Calculation Info", padding=5); ff.pack(fill="x", pady=(0, 5))
        ttk.Label(ff, text="Standard: Nearest Neighbor (NN) | [Na+]=50mM, [Primer]=0.5uM", font=("TkDefaultFont", 9, "italic"), foreground="blue").pack()
        
        # 메인 컨테이너 (왼쪽 입력 / 오른쪽 결과)
        main_split = ttk.PanedWindow(f, orient="horizontal")
        main_split.pack(fill="both", expand=True)
        
        left_f = ttk.Frame(main_split)
        main_split.add(left_f, weight=1)
        
        in_f = ttk.LabelFrame(left_f, text="Primer Registration", padding=10); in_f.pack(fill="x")
        
        r0 = ttk.Frame(in_f); r0.pack(fill="x")
        ttk.Label(r0, text="Name:").pack(side="left"); self.pn = ttk.Entry(r0, width=20); self.pn.pack(side="left", padx=5)
        
        self.pt_v = tk.StringVar(value="Fwd"); ttk.Radiobutton(r0, text="Fwd", variable=self.pt_v, value="Fwd").pack(side="left")
        ttk.Radiobutton(r0, text="Rev", variable=self.pt_v, value="Rev").pack(side="left", padx=(0, 15))
        
        self.pus_v = tk.StringVar(value="PCR"); ttk.Radiobutton(r0, text="PCR", variable=self.pus_v, value="PCR").pack(side="left")
        ttk.Radiobutton(r0, text="RT", variable=self.pus_v, value="RT-PCR").pack(side="left")

        # Sequence Row
        r1 = ttk.Frame(in_f); r1.pack(fill="x", pady=5)
        ttk.Label(r1, text="Sequence:").pack(side="left")
        self.pf = ttk.Entry(r1)
        self.pf.pack(side="left", fill="x", expand=True, padx=5)
        self.len_l = ttk.Label(r1, text="0 bp", foreground="gray")
        self.len_l.pack(side="left")
        self.pf.bind("<KeyRelease>", self.on_seq_change)

        # Binding / Action Buttons
        r2 = ttk.Frame(in_f); r2.pack(fill="x")
        ttk.Label(r2, text="Binding Length:").pack(side="left")
        self.pbl = ttk.Spinbox(r2, from_=1, to=200, width=5); self.pbl.set(20); self.pbl.pack(side="left", padx=5)
        
        ttk.Button(r2, text="Selection Tool", command=self.open_selection_tool).pack(side="left", padx=5)
        
        # Validation & Save
        self.btn_val = ttk.Button(r2, text="✅ Validate", command=self.validate_primer)
        self.btn_val.pack(side="right", padx=5)
        self.btn_save = ttk.Button(r2, text="💾 Save", command=self.add_p, state="disabled")
        self.btn_save.pack(side="right", padx=5)
        
        # Analysis Result View (Right Side)
        self.right_f = ttk.LabelFrame(main_split, text="Analysis Results", padding=10)
        main_split.add(self.right_f, weight=1)
        self.val_text = scrolledtext.ScrolledText(self.right_f, font=("Courier New", 9), state="disabled")
        self.val_text.pack(fill="both", expand=True)

        # Bottom Tree Section
        bottom_f = ttk.Frame(left_f)
        bottom_f.pack(fill="both", expand=True, pady=10)
        
        self.search = UIHelper.add_search_bar(bottom_f, self.filter_p, self.remove_dups)
        
        # Use PanedWindow to allow resizing between Fwd and Rev tables
        self.p_pane = ttk.PanedWindow(bottom_f, orient="horizontal")
        self.p_pane.pack(fill="both", expand=True, pady=5)
        
        self.fwd_t = self.make_tree(self.p_pane, "Forward", 0)
        self.rev_t = self.make_tree(self.p_pane, "Reverse", 1)
        
        # Add the LabelFrame (parent of the tree's container) to the PanedWindow
        # self.fwd_t is Treeview -> self.fwd_t.master is Frame -> self.fwd_t.master.master is LabelFrame
        self.p_pane.add(self.fwd_t.master.master, weight=1)
        self.p_pane.add(self.rev_t.master.master, weight=1)

    def on_seq_change(self, event=None):
        self.update_len_display()
        self.validated = False
        self.btn_save.config(state="disabled")

    def validate_primer(self):
        seq = self.pf.get().strip().upper()
        if not seq:
            messagebox.showwarning("Warning", "Enter sequence first.")
            return
        
        try:
            # Hetero-dimer 분석을 위해 반대쪽 프라이머 하나를 가져오거나 함 (옵션)
            # 여기서는 현재 입력된 것만 분석
            res = get_primer_analysis(seq)
            
            self.val_text.config(state="normal")
            self.val_text.delete("1.0", tk.END)
            
            # 1. Basic Properties
            self.val_text.insert(tk.END, "=== 1. Basic Properties ===\n")
            for k, v in res['basic'].items():
                self.val_text.insert(tk.END, f"{k:15}: {v}\n")
            
            # 2. Hairpin
            self.val_text.insert(tk.END, "\n=== 2. Hairpin Analysis ===\n")
            h = res['hairpin']
            color = "RED" if h['warning'] else "BLACK"
            self.val_text.insert(tk.END, f"Delta G: {h['dg']} kcal/mol\n", "warn" if h['warning'] else "")
            self.val_text.insert(tk.END, f"Tm: {h['tm']} C\n")
            if h['warning']:
                self.val_text.insert(tk.END, "⚠️ WARNING: Hairpin formation likely!\n", "warn")
                
            # 3. Homo-dimer
            self.val_text.insert(tk.END, "\n=== 3. Homo-dimer Analysis ===\n")
            hd = res['homodimer']
            self.val_text.insert(tk.END, f"Max Delta G: {hd['dg']} kcal/mol\n")
            self.val_text.insert(tk.END, f"Structure:\n{hd['structure']}\n")
            
            self.val_text.tag_config("warn", foreground="red", font=("Courier New", 9, "bold"))
            self.val_text.config(state="disabled")
            
            self.validated = True
            self.last_validated_seq = seq
            self.btn_save.config(state="normal")
            
        except Exception as e:
            messagebox.showerror("Analysis Error", f"Failed to analyze primer:\n{e}\n\nMake sure 'primer3-py' is installed.")
    def open_selection_tool(self):
        seq = self.pf.get().strip().upper()
        if not seq:
            messagebox.showinfo("Hint", "Enter sequence first.")
            return
        
        tw = tk.Toplevel(self.parent)
        tw.title("Binding Site Selector")
        tw.geometry("500x300")
        
        ttk.Label(tw, text="Drag to select the Binding Part (usually the 3' end):", padding=10).pack()
        txt = tk.Text(tw, height=5, font=("Courier New", 12), wrap="char")
        txt.pack(fill="both", expand=True, padx=10, pady=5)
        txt.insert("1.0", seq)
        
        def set_binding():
            try:
                sel = txt.tag_ranges("sel")
                if sel:
                    start_idx = len(txt.get("1.0", sel[0]))
                    end_idx = len(txt.get("1.0", sel[1]))
                    # Primer는 뒤쪽(3')이 붙는 경우가 기본이므로, 선택된 영역이 뒤에서 몇 bp인지 계산
                    # 하지만 여기서는 사용자가 선택한 그 영역 자체의 길이를 binding_len으로 사용
                    # (구조적으로 Primer 클래스는 뒤에서부터 binding_len 만큼을 binding site로 간주함)
                    selected_len = end_idx - start_idx
                    if end_idx != len(seq):
                        if not messagebox.askyesno("Warning", "Usually, the binding site is at the very 3' end. Continue?"):
                            return
                    self.pbl.set(selected_len)
                    tw.destroy()
                else:
                    messagebox.showwarning("No Selection", "Please drag to select a region.")
            except Exception as e:
                messagebox.showerror("Error", str(e))

        ttk.Button(tw, text="Set Selected as Binding Site", command=set_binding).pack(pady=10)

    def make_tree(self, parent, title, col):
        fr = ttk.LabelFrame(parent, text=title)
        # We don't use grid here because it will be added to a PanedWindow
        table_id = "primer_fwd" if title == "Forward" else "primer_rev"
        ft, t = UIHelper.create_scrolled_tree(fr, ("N", "U", "Tm_B", "Tm_T", "S", "L"), 
                                             ("Name", "Usage", "Tm(Bind)", "Tm(Total)", "Seq", "Len"),
                                             lib=self.lib, table_id=table_id)
        ft.pack(fill="both", expand=True); t.column("Tm_B", width=60); t.column("Tm_T", width=60)
        bf = ttk.Frame(fr); bf.pack(fill="x")
        ttk.Button(bf, text="Edit", command=lambda: self.start_edit(t)).pack(side="left")
        ttk.Button(bf, text="Del", command=lambda: self.del_p(t)).pack(side="right")
        return t

    def add_p(self):
        try:
            n = self.pn.get().strip()
            s = self.pf.get().strip().upper()
            if not n or not s:
                return
            binding_len = int(float(self.pbl.get()))
            t = self.pt_v.get()
            u = self.pus_v.get()
            p_obj = Primer(n, s, binding_len, t, u)
            (self.lib.fwd_primers if t == "Fwd" else self.lib.rev_primers)[n] = p_obj
            self.lib.save()
            self.refresh_cb()
            # 저장 후 입력 필드 초기화
            self.pn.delete(0, tk.END)
            self.pf.delete(0, tk.END)
            self.update_len_display()
        except Exception as e:
            from tkinter import messagebox
            messagebox.showerror("Save Error", f"Failed to save primer:\n{e}")

    def filter_p(self, q):
        for t, d in [(self.fwd_t, self.lib.fwd_primers), (self.rev_t, self.lib.rev_primers)]:
            t.delete(*t.get_children())
            for k, v in d.items():
                if not q or q.upper() in k.upper(): 
                    t.insert("", "end", values=(v.name, v.usage, v.tm, v.tm_total, v.full_sequence, len(v.full_sequence)))

    def start_edit(self, t): self.edit_widgets = UIHelper.setup_inline_edit(t, {"U":["PCR", "RT-PCR"], "S":None}, {"U":1, "S":3}, lambda mid: self.save_edit(t, mid))
    def save_edit(self, t, mid):
        vals = [w.get() for w in self.edit_widgets]; name = t.item(mid)['values'][0]; d = self.lib.fwd_primers if t==self.fwd_t else self.lib.rev_primers
        if name in d: p = d[name]; d[name] = Primer(name, vals[1], p.binding_len, p.p_type, vals[0]); self.lib.save(); self.refresh_cb()

    def remove_dups(self):
        def c(d):
            seen, dups = set(), []
            for k, v in d.items():
                if v.full_sequence in seen: dups.append(k)
                else: seen.add(v.full_sequence)
            for k in dups: del d[k]
            return len(dups)
        cnt = c(self.lib.fwd_primers) + c(self.lib.rev_primers)
        if cnt: self.lib.save(); self.refresh_cb(); messagebox.showinfo("Done", f"{cnt} dups removed")

    def upload_excel(self):
        p = filedialog.askopenfilename(filetypes=[("Excel", "*.xlsx")])
        if p:
            try:
                df = pd.read_excel(p); df.columns = [str(c).strip().lower() for c in df.columns]
                cnt = 0
                for _, r in df.iterrows():
                    n, s = str(r.get('name','')).strip(), str(r.get('sequence','')).strip()
                    if n=='nan' or not n or not s: continue
                    pt = "Fwd" if 'fwd' in str(r.get('type','')).lower() else "Rev"
                    us = "RT-PCR" if 'rt' in str(r.get('usage','pcr')).lower() else "PCR"
                    bl = int(r.get('bindinglength', 20))
                    (self.lib.fwd_primers if pt=="Fwd" else self.lib.rev_primers)[n] = Primer(n, s, bl, pt, us); cnt += 1
                self.lib.save(); self.refresh_cb()
            except Exception as e: messagebox.showerror("Error", str(e))
