import tkinter as tk
from tkinter import ttk, messagebox, filedialog, scrolledtext
import csv, os
from core.models import Species, AnnotationRef, SequenceItem
from gui.gui_components import UIHelper, FeatureEditor, SearchableCombobox, SnapGeneViewer
from core.utils import parse_annotation_file, get_primer_analysis

class PrimerValidationTab:
    def __init__(self, parent, lib, refresh_callback):
        self.parent, self.lib, self.refresh_cb = parent, lib, refresh_callback
        self.setup_ui()

    def setup_ui(self):
        f = ttk.Frame(self.parent, padding=20); f.pack(fill="both", expand=True)
        
        # 1. 상단 입력부
        in_f = ttk.LabelFrame(f, text="Primer Input", padding=10); in_f.pack(fill="x")
        
        r1 = ttk.Frame(in_f); r1.pack(fill="x", pady=2)
        ttk.Label(r1, text="Primer 1 (Fwd):", width=15).pack(side="left")
        self.p1_seq = ttk.Entry(r1); self.p1_seq.pack(side="left", fill="x", expand=True, padx=5)
        self.p1_lib = SearchableCombobox(r1, width=20); self.p1_lib.pack(side="left", padx=5)
        self.p1_lib.bind("<<ComboboxSelected>>", lambda e: self.load_from_lib(1))
        
        r2 = ttk.Frame(in_f); r2.pack(fill="x", pady=2)
        ttk.Label(r2, text="Primer 2 (Rev):", width=15).pack(side="left")
        self.p2_seq = ttk.Entry(r2); self.p2_seq.pack(side="left", fill="x", expand=True, padx=5)
        self.p2_lib = SearchableCombobox(r2, width=20); self.p2_lib.pack(side="left", padx=5)
        self.p2_lib.bind("<<ComboboxSelected>>", lambda e: self.load_from_lib(2))
        
        btn_f = ttk.Frame(in_f); btn_f.pack(fill="x", pady=5)
        ttk.Button(btn_f, text="🔍 Analyze & Validate", command=self.run_analysis).pack(side="right")
        ttk.Label(btn_f, text="* Compare two primers side-by-side", font=("TkDefaultFont", 8, "italic")).pack(side="left")

        # 2. 중간 분석 결과 (좌우 분할)
        mid_panel = ttk.PanedWindow(f, orient="horizontal")
        mid_panel.pack(fill="both", expand=True, pady=10)
        
        # Primer 1 Result
        p1_res_f = ttk.LabelFrame(mid_panel, text="Primer 1 Individual Analysis", padding=5)
        self.res_p1 = scrolledtext.ScrolledText(p1_res_f, font=("Courier New", 9), height=15)
        self.res_p1.pack(fill="both", expand=True)
        mid_panel.add(p1_res_f, weight=1)
        
        # Primer 2 Result
        p2_res_f = ttk.LabelFrame(mid_panel, text="Primer 2 Individual Analysis", padding=5)
        self.res_p2 = scrolledtext.ScrolledText(p2_res_f, font=("Courier New", 9), height=15)
        self.res_p2.pack(fill="both", expand=True)
        mid_panel.add(p2_res_f, weight=1)

        # 3. 하단 통합 결과 (Hetero-dimer & Recommendation)
        comb_res_f = ttk.LabelFrame(f, text="Pairwise Analysis & Recommendation", padding=5)
        comb_res_f.pack(fill="x", pady=(5, 0))
        self.res_comb = scrolledtext.ScrolledText(comb_res_f, font=("Courier New", 10), height=8)
        self.res_comb.pack(fill="both", expand=True)
        
        # Tag Configs for all text widgets
        for w in [self.res_p1, self.res_p2, self.res_comb]:
            w.tag_config("title", font=("Courier New", 10, "bold"), foreground="#2c3e50")
            w.tag_config("warn", font=("Courier New", 9, "bold"), foreground="red")
            w.tag_config("header", font=("Courier New", 9, "bold"), background="#e9ecef")
            w.tag_config("score", font=("Courier New", 11, "bold"), foreground="blue")

    def load_from_lib(self, target):
        combobox = self.p1_lib if target == 1 else self.p2_lib
        entry = self.p1_seq if target == 1 else self.p2_seq
        name = combobox.get()
        p_obj = self.lib.fwd_primers.get(name) or self.lib.rev_primers.get(name)
        if p_obj:
            entry.delete(0, tk.END)
            entry.insert(0, p_obj.full_sequence)

    def refresh_lists(self):
        all_p = sorted(list(self.lib.fwd_primers.keys()) + list(self.lib.rev_primers.keys()))
        self.p1_lib['values'] = all_p
        self.p2_lib['values'] = all_p

    def run_analysis(self):
        s1 = self.p1_seq.get().strip().upper()
        s2 = self.p2_seq.get().strip().upper()
        if not s1:
            messagebox.showwarning("Warning", "Please enter at least Primer 1 sequence.")
            return
            
        try:
            res = get_primer_analysis(s1, s2 if s2 else None)
            
            # Reset all views
            for w in [self.res_p1, self.res_p2, self.res_comb]:
                w.config(state="normal")
                w.delete("1.0", tk.END)

            def print_to_widget(w, data):
                w.insert(tk.END, " [1. Basic Properties]\n", "header")
                for k, v in data['basic'].items():
                    w.insert(tk.END, f"  {k:10}: {v}\n")
                
                w.insert(tk.END, "\n [2. Hairpin Analysis]\n", "header")
                h = data['hairpin']
                w.insert(tk.END, f"  Delta G: {h['dg']} kcal/mol\n")
                w.insert(tk.END, f"  Delta H: {h['dh']} kcal/mol\n")
                w.insert(tk.END, f"  Delta S: {h['ds']} kcal/K/mol\n")
                w.insert(tk.END, f"  Tm     : {h['tm']} C\n")
                if h['warning']: w.insert(tk.END, "  ⚠️ WARNING: Stable hairpin!\n", "warn")
                
                w.insert(tk.END, "\n [3. Self-Dimer Analysis]\n", "header")
                hd = data['homodimer']
                w.insert(tk.END, f"  Max dG : {hd['dg']} kcal/mol\n")
                w.insert(tk.END, f"  Delta H: {hd['dh']} kcal/mol\n")
                w.insert(tk.END, f"  Delta S: {hd['ds']} kcal/K/mol\n")
                w.insert(tk.END, f"  Structure:\n{hd['structure']}\n")

            # Update P1 and P2 side-by-side
            print_to_widget(self.res_p1, res['p1'])
            if 'p2' in res:
                print_to_widget(self.res_p2, res['p2'])
                
                # Update Combined section
                c = self.res_comb
                c.insert(tk.END, " >>> HETERO-DIMER ANALYSIS <<<\n", "header")
                htd = res['heterodimer']
                c.insert(tk.END, f" Max Delta G: {htd['dg']} kcal/mol\n")
                c.insert(tk.END, f" Delta H     : {htd['dh']} kcal/mol\n")
                c.insert(tk.END, f" Delta S     : {htd['ds']} kcal/K/mol\n")
                c.insert(tk.END, f" Structure:\n{htd['structure']}\n")
                
                # Recommendation Score
                score = 100
                reasons = []
                min_dg = min(res['p1']['hairpin']['dg'], res['p1']['homodimer']['dg'],
                             res['p2']['hairpin']['dg'], res['p2']['homodimer']['dg'],
                             res['heterodimer']['dg'])
                if min_dg < -9.0: score -= 50; reasons.append("Critical secondary structure")
                elif min_dg < -5.0: score -= 20; reasons.append("Moderate secondary structure")
                
                tm_diff = abs(res['p1']['basic']['Tm'] - res['p2']['basic']['Tm'])
                if tm_diff > 5.0: score -= 30; reasons.append(f"High Tm difference ({tm_diff}C)")
                elif tm_diff > 2.0: score -= 10; reasons.append(f"Minor Tm difference ({tm_diff}C)")

                c.insert(tk.END, "-"*60 + "\n")
                c.insert(tk.END, f" FINAL SCORE: {max(0, score)} / 100  |  ", "score")
                
                if score >= 90: grade = "⭐⭐⭐⭐⭐ (Excellent)"
                elif score >= 70: grade = "⭐⭐⭐⭐ (Good)"
                elif score >= 40: grade = "⭐⭐ (Cautious)"
                else: grade = "⭐ (Not Recommended)"
                
                c.insert(tk.END, f"Grade: {grade}\n", "score")
                if reasons: c.insert(tk.END, " Issues: " + ", ".join(reasons) + "\n", "warn")
            else:
                self.res_p2.insert(tk.END, "\n\n  (No Primer 2 entered)")
                self.res_comb.insert(tk.END, " Enter Primer 2 for Pairwise Analysis & Recommendation.")

            for w in [self.res_p1, self.res_p2, self.res_comb]:
                w.config(state="disabled")
        except Exception as e:
            messagebox.showerror("Error", str(e))
import os

class SpeciesTab:
    def __init__(self, parent, lib, refresh_callback):
        self.parent, self.lib, self.refresh_cb = parent, lib, refresh_callback
        self.setup_ui()

    def setup_ui(self):
        f = ttk.Frame(self.parent, padding=10); f.pack(fill="both", expand=True)
        in_f = ttk.LabelFrame(f, text="Species Manager", padding=10); in_f.pack(fill="x")
        ttk.Label(in_f, text="Name:").pack(side="left"); self.sn = ttk.Entry(in_f); self.sn.pack(side="left", padx=5)
        ttk.Label(in_f, text="Desc:").pack(side="left"); self.sd = ttk.Entry(in_f); self.sd.pack(side="left", padx=5)
        ttk.Button(in_f, text="Save", command=self.add_s).pack(side="left")
        
        self.search = UIHelper.add_search_bar(f, self.filter, lambda: None)
        fr_t, self.tree = UIHelper.create_scrolled_tree(f, ("N", "D"), ("Name", "Description"), lib=self.lib, table_id="species_main")
        fr_t.pack(fill="both", expand=True)
        
        bf = ttk.Frame(f); bf.pack(fill="x")
        ttk.Button(bf, text="Edit", command=self.edit_load).pack(side="left")
        ttk.Button(bf, text="Delete", command=self.del_s).pack(side="right")

    def edit_load(self):
        sel = self.tree.selection()
        if sel:
            s = self.lib.species[self.tree.item(sel[0])['values'][0]]
            self.sn.delete(0, tk.END); self.sn.insert(0, s.name); self.sd.delete(0, tk.END); self.sd.insert(0, s.description)

    def add_s(self):
        n = self.sn.get().strip()
        if n: self.lib.species[n] = Species(n, self.sd.get()); self.lib.save(); self.refresh_cb()

    def del_s(self):
        for s in self.tree.selection(): del self.lib.species[self.tree.item(s)['values'][0]]
        self.lib.save(); self.refresh_cb()

    def filter(self, q):
        self.tree.delete(*self.tree.get_children())
        for k, v in self.lib.species.items():
            if not q or q.upper() in k.upper(): self.tree.insert("", "end", values=(v.name, v.description))

from core.providers import ProviderRegistry

class AnnotationRefTab:
    def __init__(self, parent, lib, refresh_callback):
        self.parent, self.lib, self.refresh_cb = parent, lib, refresh_callback
        self.setup_ui()

    def setup_ui(self):
        f = ttk.Frame(self.parent, padding=10); f.pack(fill="both", expand=True)
        in_f = ttk.LabelFrame(f, text="Annotation Reference Manager", padding=10); in_f.pack(fill="x")
        ttk.Button(in_f, text="📁 Load GFF3/GBK", command=self.load_file).pack(side="left")
        ttk.Button(in_f, text="?", width=3, command=self.show_help).pack(side="left", padx=5)
        
        self.search = UIHelper.add_search_bar(f, self.filter, lambda: None)
        fr_t, self.tree = UIHelper.create_scrolled_tree(f, ("N", "F", "Sp", "S"), ("Name", "Features", "Species", "Source"), lib=self.lib, table_id="ann_ref_main")
        fr_t.pack(fill="both", expand=True)

        # Add Species/Provider edit area
        edit_f = ttk.Frame(f, padding=5); edit_f.pack(fill="x")
        ttk.Label(edit_f, text="Species:").pack(side="left")
        self.sp_c = ttk.Combobox(edit_f, values=ProviderRegistry.list_species())
        self.sp_c.pack(side="left", padx=5)
        ttk.Button(edit_f, text="Update Selected Species", command=self.update_species).pack(side="left")

        ttk.Button(f, text="Delete", command=self.del_r).pack(pady=5)

    def show_help(self): messagebox.showinfo("Help", "Use .gff3 or .gbk files. IDs must match FASTA sequence IDs.")
    
    def load_file(self):
        p = filedialog.askopenfilename(filetypes=[("Annotation", "*.gff *.gff3 *.gbk *.gb")])
        if p:
            feats = parse_annotation_file(p); n = os.path.basename(p)
            # Default to Yeast if SGD mentioned, else Generic
            species = "Saccharomyces cerevisiae" if "sgd" in p.lower() or "yeast" in p.lower() else "Generic"
            provider_obj = ProviderRegistry.get_by_species(species)
            self.lib.ann_refs[n] = AnnotationRef(n, feats, p, species, provider_obj.name)
            self.lib.save(); self.refresh_cb()

    def update_species(self):
        sel = self.tree.selection()
        if not sel: return
        name = self.tree.item(sel[0])['values'][0]
        if name in self.lib.ann_refs:
            sp = self.sp_c.get()
            self.lib.ann_refs[name].species = sp
            self.lib.ann_refs[name].provider = ProviderRegistry.get_by_species(sp).name
            self.lib.save()
            self.refresh_cb()

    def del_r(self):
        for s in self.tree.selection(): 
            val = self.tree.item(s)['values']
            if val and val[0] in self.lib.ann_refs:
                del self.lib.ann_refs[val[0]]
        self.lib.save(); self.refresh_cb()

    def filter(self, q):
        self.tree.delete(*self.tree.get_children())
        for k, v in self.lib.ann_refs.items():
            if not q or q.upper() in k.upper(): 
                self.tree.insert("", "end", values=(v.name, len(v.features), getattr(v, 'species', 'Generic'), v.source_file))

class SeqViewTab:
    def __init__(self, parent, lib, refresh_callback, data_dict, category):
        self.parent, self.lib, self.refresh_cb, self.data, self.cat = parent, lib, refresh_callback, data_dict, category
        self.setup_ui()

    def setup_ui(self):
        f = ttk.Frame(self.parent, padding=10); f.pack(fill="both", expand=True)
        bf = ttk.LabelFrame(f, text="Update Metadata", padding=10); bf.pack(fill="x")
        self.ae_n = ttk.Entry(bf, width=30); self.ae_n.pack(side="left", padx=5)
        self.ae_m = SearchableCombobox(bf); self.ae_m.pack(side="left", padx=5)
        self.ae_t = SearchableCombobox(bf); self.ae_t.pack(side="left", padx=5)
        ttk.Button(bf, text="Update Info", command=self.upd).pack(side="left")

        # Export buttons
        exp_f = ttk.Frame(f); exp_f.pack(fill="x", pady=2)
        ttk.Button(exp_f, text="Export FASTA", command=self.export_fasta).pack(side="left", padx=2)
        ttk.Button(exp_f, text="Export CSV", command=self.export_csv).pack(side="left", padx=2)

        self.search = UIHelper.add_search_bar(f, self.filter, lambda: self.remove_dups())
        if self.cat == "amplicon":
            cols, headings = ("N", "L", "Tmp", "M", "T"), ("Name", "Len", "Template", "Marker", "Tag")
        elif self.cat == "digest":
            cols, headings = ("N", "L", "Src", "Enz"), ("Name", "Len", "Source", "Enzymes")
        else:
            cols, headings = ("N", "L", "M", "T"), ("Name", "Len", "Marker", "Tag")
        fr_t, self.tree = UIHelper.create_scrolled_tree(f, cols, headings, lib=self.lib, table_id=f"seq_view_{self.cat}")
        fr_t.pack(fill="x"); self.tree.bind("<<TreeviewSelect>>", lambda e: self.on_sel())

        self.det = SnapGeneViewer(f, lib_manager=self.lib)
        self.det.pack(fill="both", expand=True, pady=5)

        self.feat_ed = FeatureEditor(f, "Features", text_widget=self.det); self.feat_ed.pack(fill="x")
        UIHelper.create_legend(f).pack(side="bottom", fill="x", pady=5)
        ttk.Button(f, text="Delete Selected", command=self.del_s).pack(side="right")

    def on_sel(self):
        sel = self.tree.selection()
        if sel:
            name = self.tree.item(sel[0])['values'][0]; item = self.data.get(name)
            if item:
                self.ae_n.delete(0, tk.END); self.ae_n.insert(0, item.name)
                self.ae_m.set(item.marker); self.ae_t.set(item.tag)
                UIHelper.render_gene_viewer(self.det, item.sequence, item.features, lib_manager=self.lib)
                self.feat_ed.set_features(item.features)

    def upd(self):
        sel = self.tree.selection()
        if sel:
            old = self.tree.item(sel[0])['values'][0]; new = self.ae_n.get().strip(); item = self.data.get(old)
            if item and new:
                item.name, item.marker, item.tag, item.features = new, self.ae_m.get(), self.ae_t.get(), self.feat_ed.get_features()
                if old != new: self.data[new] = item; del self.data[old]
                self.lib.save(); self.refresh_cb()

    def del_s(self):
        for s in self.tree.selection():
            del self.data[self.tree.item(s)['values'][0]]
        self.lib.save()
        self.refresh_cb()
        self.det.show_message("")

    def export_fasta(self):
        path = filedialog.asksaveasfilename(
            defaultextension=".fasta",
            filetypes=[("FASTA", "*.fasta *.fa"), ("All", "*.*")],
            title="Export as FASTA")
        if not path:
            return
        sel = self.tree.selection()
        items = ([self.data[self.tree.item(s)['values'][0]] for s in sel]
                 if sel else list(self.data.values()))
        with open(path, 'w', encoding='utf-8') as fh:
            for item in items:
                fh.write(f">{item.name}\n")
                seq = item.sequence
                for i in range(0, len(seq), 60):
                    fh.write(seq[i:i+60] + "\n")
        messagebox.showinfo("Export", f"{len(items)}개 서열을 저장했습니다.\n{path}")

    def export_csv(self):
        path = filedialog.asksaveasfilename(
            defaultextension=".csv",
            filetypes=[("CSV", "*.csv"), ("All", "*.*")],
            title="Export as CSV")
        if not path:
            return
        sel = self.tree.selection()
        items = ([self.data[self.tree.item(s)['values'][0]] for s in sel]
                 if sel else list(self.data.values()))
        with open(path, 'w', newline='', encoding='utf-8-sig') as fh:
            writer = csv.writer(fh)
            writer.writerow(["Name", "Length", "Topology", "Marker", "Tag", "Template", "Sequence"])
            for item in items:
                writer.writerow([
                    item.name, len(item.sequence),
                    getattr(item, 'topology', ''),
                    getattr(item, 'marker', ''),
                    getattr(item, 'tag', ''),
                    getattr(item, 'template_name', ''),
                    item.sequence])
        messagebox.showinfo("Export", f"{len(items)}개 서열을 저장했습니다.\n{path}")

    def filter(self, q):
        self.tree.delete(*self.tree.get_children())
        for k, v in self.data.items():
            if not q or q.upper() in k.upper():
                if self.cat == "amplicon":
                    row = (v.name, f"{len(v.sequence)}bp", v.template_name, v.marker, v.tag)
                elif self.cat == "digest":
                    row = (v.name, f"{len(v.sequence)}bp", v.template_name, getattr(v, 'enzymes', ''))
                else:
                    row = (v.name, f"{len(v.sequence)}bp", v.marker, v.tag)
                self.tree.insert("", "end", values=row)

    def remove_dups(self):
        seen, dups = set(), []
        for k, v in self.data.items():
            if v.sequence in seen: dups.append(k)
            else: seen.add(v.sequence)
        for k in dups: del self.data[k]
        if dups: self.lib.save(); self.refresh_cb()
