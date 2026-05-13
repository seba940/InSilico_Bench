import tkinter as tk
from tkinter import ttk, messagebox
import traceback
from core.models import SequenceItem
from gui.gui_components import UIHelper, SearchableCombobox, SnapGeneViewer
from core.utils import (get_rc, simulate_hr_logic_with_features,
                        calculate_extension_time, smart_find_binding,
                        adjust_features, POLYMERASE_RATES)


class PCRTab:
    def __init__(self, parent, lib, refresh_callback):
        self.parent, self.lib, self.refresh_cb = parent, lib, refresh_callback
        self.setup_ui()

    def setup_ui(self):
        f = ttk.Frame(self.parent, padding=20); f.pack(fill="both", expand=True)
        sf = ttk.Frame(f); sf.pack(fill="x")
        ttk.Label(sf, text="Template:").grid(row=0, column=0, sticky="w")
        self.tmp = SearchableCombobox(sf); self.tmp.grid(row=0, column=1, sticky="ew", pady=2)
        ttk.Label(sf, text="Fwd Primer:").grid(row=1, column=0, sticky="w")
        self.fwd = SearchableCombobox(sf); self.fwd.grid(row=1, column=1, sticky="ew", pady=2)
        ttk.Label(sf, text="Rev Primer:").grid(row=2, column=0, sticky="w")
        self.rev = SearchableCombobox(sf); self.rev.grid(row=2, column=1, sticky="ew", pady=2)
        ttk.Label(sf, text="Polymerase:").grid(row=3, column=0, sticky="w")
        self.poly = ttk.Combobox(sf, values=list(POLYMERASE_RATES.keys()), state="readonly", width=12)
        self.poly.set("Taq"); self.poly.grid(row=3, column=1, sticky="w", pady=2)
        self.poly.bind("<<ComboboxSelected>>", lambda e: self.upd_info())
        sf.columnconfigure(1, weight=1)

        if0 = ttk.LabelFrame(f, text="PCR Design Info", padding=10); if0.pack(fill="x", pady=10)
        self.tm_l = ttk.Label(if0, text="Annealing Temp (Ta): - ℃", font=("TkDefaultFont", 10, "bold")); self.tm_l.pack(side="left", padx=20)
        self.ex_l = ttk.Label(if0, text="Extension Time: -", font=("TkDefaultFont", 10)); self.ex_l.pack(side="left", padx=20)
        self.ln_l = ttk.Label(if0, text="Product Length: - bp", font=("TkDefaultFont", 10, "bold"), foreground="blue"); self.ln_l.pack(side="left", padx=20)

        self.tmp.bind("<<ComboboxSelected>>", lambda e: self.upd_info())
        self.fwd.bind("<<ComboboxSelected>>", lambda e: self.upd_info())
        self.rev.bind("<<ComboboxSelected>>", lambda e: self.upd_info())
        ttk.Button(f, text="⚡ Run PCR Simulation", command=self.run).pack(pady=5)

        self.res = SnapGeneViewer(f, lib_manager=self.lib)
        self.res.pack(fill="both", expand=True)

        ttk.Button(f, text="Save Amplicon to Library", command=self.save).pack(pady=5)
        UIHelper.create_legend(f).pack(side="bottom", fill="x", pady=5)

    def upd_info(self):
        f = self.lib.fwd_primers.get(self.fwd.get())
        r = self.lib.rev_primers.get(self.rev.get())
        if f and r:
            ta = min(f.tm, r.tm)
            self.tm_l.config(text=f"Annealing Temp (Ta): {ta} ℃")

            t_sel = self.tmp.get().split("] ")[-1] if "]" in self.tmp.get() else self.tmp.get()
            template_obj = self.lib.templates.get(t_sel) or self.lib.digests.get(t_sel) or self.lib.recombinants.get(t_sel)
            if template_obj:
                UIHelper.render_gene_viewer(self.res, template_obj.sequence, template_obj.features,
                                            primers_to_show=[f, r],
                                            lib_manager=self.lib)

    def run(self):
        t_sel = self.tmp.get().split("] ")[-1] if "]" in self.tmp.get() else self.tmp.get()
        template_obj = self.lib.templates.get(t_sel) or self.lib.digests.get(t_sel) or self.lib.recombinants.get(t_sel)
        f_primer = self.lib.fwd_primers.get(self.fwd.get())
        r_primer = self.lib.rev_primers.get(self.rev.get())

        if not (template_obj and f_primer and r_primer):
            messagebox.showwarning("Warning", "Select Template and both Primers.")
            return

        try:
            idx_f, len_f = smart_find_binding(template_obj.sequence, f_primer.binding, "Fwd", template_obj.topology)
            idx_r, len_r = smart_find_binding(template_obj.sequence, r_primer.binding, "Rev", template_obj.topology)

            if idx_f == -1 or idx_r == -1:
                messagebox.showerror("Failed", "One or both primers failed to bind to the template.")
                return

            if template_obj.topology == "Circular":
                if idx_r < idx_f:
                    mid_seq = template_obj.sequence[idx_f + len_f:] + template_obj.sequence[:idx_r]
                else:
                    mid_seq = template_obj.sequence[idx_f + len_f:idx_r]
            else:
                if idx_r < idx_f:
                    msg = f"Primers are facing wrong directions for Linear PCR.\n\nFwd binding: {idx_f} ~ {idx_f+len_f}\nRev binding: {idx_r} ~ {idx_r+len_r}"
                    idx_f_alt, _ = smart_find_binding(template_obj.sequence, r_primer.binding, "Fwd", template_obj.topology)
                    idx_r_alt, _ = smart_find_binding(template_obj.sequence, f_primer.binding, "Rev", template_obj.topology)
                    if idx_f_alt != -1 and idx_r_alt != -1 and idx_f_alt < idx_r_alt:
                        msg += "\n\nTip: It seems swapping Fwd and Rev primers might work."
                    messagebox.showerror("Failed", msg)
                    return
                mid_seq = template_obj.sequence[idx_f + len_f:idx_r]

            amp_seq = f_primer.full_sequence + mid_seq + get_rc(r_primer.full_sequence)

            t_len = len(template_obj.sequence)
            end_idx = idx_r + len_r
            if template_obj.topology == "Circular":
                end_idx = end_idx % t_len
            new_feats = adjust_features(template_obj.features, len(f_primer.overhang), idx_f, end_idx, template_obj.topology, t_len)
            new_feats.append({"label": f_primer.name, "start": 0, "end": len(f_primer.full_sequence), "type": "Homology Arm"})
            new_feats.append({"label": r_primer.name, "start": len(amp_seq)-len(r_primer.full_sequence), "end": len(amp_seq), "type": "Homology Arm"})

            self.last_amp = SequenceItem(
                name=f"{f_primer.name}-{r_primer.name}",
                sequence=amp_seq,
                category="amplicon",
                features=new_feats,
                template_name=t_sel
            )

            UIHelper.render_gene_viewer(self.res, amp_seq, new_feats,
                                        primers_to_show=[f_primer, r_primer],
                                        lib_manager=self.lib)

            self.ln_l.config(text=f"Product Length: {len(amp_seq)} bp")
            poly = self.poly.get()
            ext_time = calculate_extension_time(len(amp_seq), poly)
            self.ex_l.config(text=f"Extension Time ({poly}): {ext_time//60}m {ext_time%60}s")

        except Exception as e:
            messagebox.showerror("Error", f"PCR Simulation failed: {str(e)}\n{traceback.format_exc()}")

    def save(self):
        if hasattr(self, 'last_amp'):
            self.lib.push_undo()
            self.lib.amplicons[self.last_amp.name] = self.last_amp
            self.lib.save()
            self.refresh_cb()
            messagebox.showinfo("Success", f"Amplicon '{self.last_amp.name}' saved to library.")


class HRTab:
    def __init__(self, parent, lib, refresh_callback):
        self.parent, self.lib, self.refresh_cb = parent, lib, refresh_callback
        self.setup_ui()

    def setup_ui(self):
        f = ttk.Frame(self.parent, padding=20); f.pack(fill="both", expand=True)
        ttk.Label(f, text="Target Genome (Kind=Genome only):").pack(anchor="w")
        self.target = SearchableCombobox(f); self.target.pack(fill="x", pady=2)
        ttk.Label(f, text="Insert (Amplicon / Digest fragment):").pack(anchor="w")
        self.amp = SearchableCombobox(f); self.amp.pack(fill="x", pady=2)

        mm_f = ttk.Frame(f); mm_f.pack(fill="x", pady=2)
        ttk.Label(mm_f, text="Max mismatches in homology arm:").pack(side="left")
        self.mm_sp = ttk.Spinbox(mm_f, from_=0, to=5, width=4)
        self.mm_sp.set(2); self.mm_sp.pack(side="left", padx=4)

        ttk.Button(f, text="⚡ Run HR Simulation", command=self.run).pack(pady=5)

        self.res = SnapGeneViewer(f, lib_manager=self.lib)
        self.res.pack(fill="both", expand=True)

        # Colony PCR verification section
        col_f = ttk.LabelFrame(f, text="Colony PCR Verification", padding=8)
        col_f.pack(fill="x", pady=5)
        col_row = ttk.Frame(col_f); col_row.pack(fill="x")
        ttk.Label(col_row, text="Check Fwd:").pack(side="left")
        self.col_fwd = SearchableCombobox(col_row); self.col_fwd.pack(side="left", fill="x", expand=True, padx=4)
        ttk.Label(col_row, text="Check Rev:").pack(side="left")
        self.col_rev = SearchableCombobox(col_row); self.col_rev.pack(side="left", fill="x", expand=True, padx=4)
        ttk.Button(col_row, text="Predict Band", command=self.verify_colony_pcr).pack(side="left", padx=4)
        self.col_result = ttk.Label(col_f, text="", font=("TkDefaultFont", 9, "bold"), foreground="darkblue")
        self.col_result.pack(anchor="w", pady=2)

        save_f = ttk.Frame(f); save_f.pack(fill="x", pady=4)
        ttk.Label(save_f, text="Save name:").pack(side="left")
        self.save_name_e = ttk.Entry(save_f, width=40)
        self.save_name_e.pack(side="left", padx=4)
        ttk.Button(save_f, text="Save Result to Library", command=self.save).pack(side="left")

        UIHelper.create_legend(f).pack(side="bottom", fill="x", pady=5)

    def _refresh_colony_primers(self):
        all_p_fwd = sorted(self.lib.fwd_primers.keys())
        all_p_rev = sorted(self.lib.rev_primers.keys())
        self.col_fwd['values'] = all_p_fwd
        self.col_rev['values'] = all_p_rev

    def run(self):
        t_sel = self.target.get().split("] ")[-1] if "]" in self.target.get() else self.target.get()
        target_obj = self.lib.templates.get(t_sel)
        amp_sel = self.amp.get()
        if amp_sel.startswith("[Digest] "):
            amp_obj = self.lib.digests.get(amp_sel[len("[Digest] "):])
        elif amp_sel.startswith("[Ligation] "):
            amp_obj = self.lib.ligations.get(amp_sel[len("[Ligation] "):])
        else:
            amp_obj = self.lib.amplicons.get(amp_sel) or self.lib.digests.get(amp_sel)

        if not (target_obj and amp_obj):
            messagebox.showwarning("Warning", "Select Target Genome and Insert.")
            return

        try:
            mm = int(self.mm_sp.get())
        except ValueError:
            mm = 2

        if target_obj and amp_obj:
            final_seq, feats, s, e = simulate_hr_logic_with_features(
                target_obj, amp_obj, max_mismatches=mm)
            if final_seq:
                self.last_hr = SequenceItem(
                    name="HR_Result",
                    sequence=final_seq,
                    category="recombinant",
                    topology=target_obj.topology,
                    features=feats,
                    marker=amp_obj.marker,
                    tag=amp_obj.tag
                )
                UIHelper.render_gene_viewer(self.res, final_seq, feats, lib_manager=self.lib)
                default_name = f"HR_{t_sel[:15]}_{amp_sel[:15]}"
                self.save_name_e.delete(0, tk.END)
                self.save_name_e.insert(0, default_name)
                self._refresh_colony_primers()
            else:
                messagebox.showerror("Error", "No homologous recombination sites found (min 10bp match at ends).")

    def verify_colony_pcr(self):
        if not hasattr(self, 'last_hr'):
            messagebox.showwarning("Warning", "먼저 HR Simulation을 실행하세요.")
            return
        fwd_name = self.col_fwd.get()
        rev_name = self.col_rev.get()
        f_primer = self.lib.fwd_primers.get(fwd_name)
        r_primer = self.lib.rev_primers.get(rev_name)
        if not (f_primer and r_primer):
            messagebox.showwarning("Warning", "Check Fwd / Rev primer를 선택하세요.")
            return
        seq = self.last_hr.sequence
        topo = self.last_hr.topology
        idx_f, len_f = smart_find_binding(seq, f_primer.binding, "Fwd", topo)
        idx_r, len_r = smart_find_binding(seq, r_primer.binding, "Rev", topo)
        if idx_f == -1 or idx_r == -1:
            self.col_result.config(
                text="Primers not found in HR result. Check primer sequences.",
                foreground="red")
            return
        if idx_r <= idx_f:
            self.col_result.config(
                text=f"Primers face wrong direction (Fwd@{idx_f}, Rev@{idx_r}).",
                foreground="red")
            return
        band = idx_r + len_r - idx_f
        self.col_result.config(
            text=f"Predicted band: {band} bp  (Fwd@{idx_f+1}, Rev@{idx_r+1})",
            foreground="darkblue")

    def save(self):
        if hasattr(self, 'last_hr'):
            n = self.save_name_e.get().strip() or "HR_Result"
            original = n; ctr = 2
            while n in self.lib.recombinants:
                n = f"{original} ({ctr})"; ctr += 1
            self.last_hr.name = n
            self.lib.push_undo()
            self.lib.recombinants[n] = self.last_hr
            self.lib.save()
            self.refresh_cb()
            messagebox.showinfo("Success", f"Recombinant '{n}' saved to library.")


class BatchPCRTab:
    """Run multiple (template, fwd, rev) combinations at once."""

    def __init__(self, parent, lib, refresh_callback):
        self.parent, self.lib, self.refresh_cb = parent, lib, refresh_callback
        self._rows = []  # list of (tmp_var, fwd_var, rev_var)
        self.setup_ui()

    def setup_ui(self):
        f = ttk.Frame(self.parent, padding=10)
        f.pack(fill="both", expand=True)

        top = ttk.Frame(f); top.pack(fill="x", pady=4)
        ttk.Label(top, text="Polymerase:").pack(side="left")
        self.poly = ttk.Combobox(top, values=list(POLYMERASE_RATES.keys()), state="readonly", width=10)
        self.poly.set("Taq"); self.poly.pack(side="left", padx=6)
        ttk.Button(top, text="+ Add Row", command=self.add_row).pack(side="left", padx=4)
        ttk.Button(top, text="Clear All", command=self.clear_rows).pack(side="left", padx=4)
        ttk.Button(top, text="⚡ Run All", command=self.run_all).pack(side="right", padx=4)
        ttk.Button(top, text="Save All Amplicons", command=self.save_all).pack(side="right", padx=4)

        # Row area (scrollable)
        row_outer = ttk.LabelFrame(f, text="PCR Reactions", padding=4)
        row_outer.pack(fill="x", pady=4)
        canvas = tk.Canvas(row_outer, height=160, highlightthickness=0)
        sb = ttk.Scrollbar(row_outer, orient="vertical", command=canvas.yview)
        canvas.configure(yscrollcommand=sb.set)
        sb.pack(side="right", fill="y")
        canvas.pack(side="left", fill="both", expand=True)
        self._row_frame = ttk.Frame(canvas)
        self._row_frame_id = canvas.create_window((0, 0), window=self._row_frame, anchor="nw")
        self._row_frame.bind("<Configure>", lambda e: canvas.configure(
            scrollregion=canvas.bbox("all")))
        canvas.bind("<Configure>", lambda e: canvas.itemconfig(
            self._row_frame_id, width=e.width))
        self._canvas = canvas

        # Header
        hdr = ttk.Frame(self._row_frame)
        hdr.pack(fill="x")
        for txt, w in [("Template", 200), ("Fwd Primer", 160), ("Rev Primer", 160), ("", 60)]:
            ttk.Label(hdr, text=txt, font=("TkDefaultFont", 8, "bold"), width=w//8).pack(side="left", padx=2)

        # Add 3 default rows
        for _ in range(3):
            self.add_row()

        # Result table
        res_f = ttk.LabelFrame(f, text="Results", padding=4)
        res_f.pack(fill="both", expand=True, pady=4)
        cols = ("template", "fwd", "rev", "length", "ta", "ext_time", "status")
        heads = ("Template", "Fwd", "Rev", "Length (bp)", "Ta (℃)", "Ext.Time", "Status")
        fr_t, self.result_tree = UIHelper.create_scrolled_tree(res_f, cols, heads)
        fr_t.pack(fill="both", expand=True)
        for col, w in zip(cols, [150, 120, 120, 90, 70, 90, 120]):
            self.result_tree.column(col, width=w)
        self._last_amplicons = []

    def _get_all_templates(self):
        tl = [f"[{v.topology}] {k}" for k, v in self.lib.templates.items()]
        dl = [f"[Digest:{v.topology}] {k}" for k, v in self.lib.digests.items()]
        rl = [f"[HR:{v.topology}] {k}" for k, v in self.lib.recombinants.items()]
        return tl + dl + rl

    def add_row(self):
        row_f = ttk.Frame(self._row_frame)
        row_f.pack(fill="x", pady=1)

        tmp_cb = SearchableCombobox(row_f)
        tmp_cb.configure(values=self._get_all_templates())
        tmp_cb.pack(side="left", fill="x", expand=True, padx=2)

        fwd_cb = SearchableCombobox(row_f)
        fwd_cb.configure(values=sorted(self.lib.fwd_primers.keys()))
        fwd_cb.pack(side="left", padx=2)

        rev_cb = SearchableCombobox(row_f)
        rev_cb.configure(values=sorted(self.lib.rev_primers.keys()))
        rev_cb.pack(side="left", padx=2)

        def _del(rf=row_f, entry=(tmp_cb, fwd_cb, rev_cb)):
            self._rows = [r for r in self._rows if r is not entry]
            rf.destroy()

        ttk.Button(row_f, text="✕", width=3, command=_del).pack(side="left", padx=2)
        entry = (tmp_cb, fwd_cb, rev_cb)
        self._rows.append(entry)

    def clear_rows(self):
        for widget in self._row_frame.winfo_children():
            if isinstance(widget, ttk.Frame) and widget != self._row_frame.winfo_children()[0]:
                widget.destroy()
        self._rows.clear()

    def refresh_combos(self):
        tmpl_vals = self._get_all_templates()
        fwd_vals = sorted(self.lib.fwd_primers.keys())
        rev_vals = sorted(self.lib.rev_primers.keys())
        for tmp_cb, fwd_cb, rev_cb in self._rows:
            tmp_cb.configure(values=tmpl_vals)
            fwd_cb.configure(values=fwd_vals)
            rev_cb.configure(values=rev_vals)

    def run_all(self):
        self.result_tree.delete(*self.result_tree.get_children())
        self._last_amplicons = []
        poly = self.poly.get()

        for tmp_cb, fwd_cb, rev_cb in self._rows:
            t_raw = tmp_cb.get().strip()
            fwd_name = fwd_cb.get().strip()
            rev_name = rev_cb.get().strip()
            if not (t_raw and fwd_name and rev_name):
                continue

            t_sel = t_raw.split("] ")[-1] if "]" in t_raw else t_raw
            template_obj = (self.lib.templates.get(t_sel)
                            or self.lib.digests.get(t_sel)
                            or self.lib.recombinants.get(t_sel))
            f_primer = self.lib.fwd_primers.get(fwd_name)
            r_primer = self.lib.rev_primers.get(rev_name)

            if not (template_obj and f_primer and r_primer):
                self.result_tree.insert("", "end", values=(
                    t_sel, fwd_name, rev_name, "-", "-", "-", "NOT FOUND"))
                continue

            try:
                idx_f, len_f = smart_find_binding(template_obj.sequence, f_primer.binding, "Fwd", template_obj.topology)
                idx_r, len_r = smart_find_binding(template_obj.sequence, r_primer.binding, "Rev", template_obj.topology)

                if idx_f == -1 or idx_r == -1 or (template_obj.topology != "Circular" and idx_r <= idx_f):
                    self.result_tree.insert("", "end", values=(
                        t_sel, fwd_name, rev_name, "-", "-", "-", "NO BINDING"))
                    continue

                if template_obj.topology == "Circular" and idx_r < idx_f:
                    mid_seq = template_obj.sequence[idx_f + len_f:] + template_obj.sequence[:idx_r]
                else:
                    mid_seq = template_obj.sequence[idx_f + len_f:idx_r]

                amp_seq = f_primer.full_sequence + mid_seq + get_rc(r_primer.full_sequence)
                ta = min(f_primer.tm, r_primer.tm)
                ext = calculate_extension_time(len(amp_seq), poly)
                ext_str = f"{ext//60}m{ext%60:02d}s"

                t_len = len(template_obj.sequence)
                end_idx = (idx_r + len_r) % t_len if template_obj.topology == "Circular" else idx_r + len_r
                new_feats = adjust_features(template_obj.features, len(f_primer.overhang), idx_f, end_idx, template_obj.topology, t_len)
                new_feats.append({"label": f_primer.name, "start": 0, "end": len(f_primer.full_sequence), "type": "Homology Arm"})
                new_feats.append({"label": r_primer.name, "start": len(amp_seq) - len(r_primer.full_sequence), "end": len(amp_seq), "type": "Homology Arm"})

                amp = SequenceItem(
                    name=f"{fwd_name}-{rev_name}",
                    sequence=amp_seq,
                    category="amplicon",
                    features=new_feats,
                    template_name=t_sel)
                self._last_amplicons.append(amp)

                self.result_tree.insert("", "end", values=(
                    t_sel, fwd_name, rev_name, len(amp_seq), ta, ext_str, "OK"))

            except Exception as ex:
                self.result_tree.insert("", "end", values=(
                    t_sel, fwd_name, rev_name, "-", "-", "-", f"ERROR: {ex}"))

    def save_all(self):
        if not self._last_amplicons:
            messagebox.showwarning("Warning", "Run All을 먼저 실행하세요.")
            return
        self.lib.push_undo()
        saved = 0
        for amp in self._last_amplicons:
            name = amp.name; ctr = 2
            while name in self.lib.amplicons:
                name = f"{amp.name} ({ctr})"; ctr += 1
            amp.name = name
            self.lib.amplicons[name] = amp
            saved += 1
        self.lib.save()
        self.refresh_cb()
        messagebox.showinfo("Saved", f"{saved}개 Amplicon이 라이브러리에 저장되었습니다.")
