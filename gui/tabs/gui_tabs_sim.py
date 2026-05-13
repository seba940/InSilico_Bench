import tkinter as tk
from tkinter import ttk, messagebox
import traceback
from core.models import SequenceItem
from gui.gui_components import UIHelper, SearchableCombobox, SnapGeneViewer
from core.utils import get_rc, simulate_hr_logic_with_features, calculate_extension_time, smart_find_binding, adjust_features


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
        self.rev = SearchableCombobox(sf); self.rev.grid(row=2, column=1, sticky="ew", pady=2); sf.columnconfigure(1, weight=1)

        if0 = ttk.LabelFrame(f, text="PCR Design Info", padding=10); if0.pack(fill="x", pady=10)
        self.tm_l = ttk.Label(if0, text="Annealing Temp (Ta): - ℃", font=("TkDefaultFont", 10, "bold")); self.tm_l.pack(side="left", padx=20)
        self.ex_l = ttk.Label(if0, text="Extension Time: - (1min/1KB)", font=("TkDefaultFont", 10)); self.ex_l.pack(side="left", padx=20)
        self.ln_l = ttk.Label(if0, text="Product Length: - bp", font=("TkDefaultFont", 10, "bold"), foreground="blue"); self.ln_l.pack(side="left", padx=20)

        self.tmp.bind("<<ComboboxSelected>>", lambda e: self.upd_info())
        self.fwd.bind("<<ComboboxSelected>>", lambda e: self.upd_info()); self.rev.bind("<<ComboboxSelected>>", lambda e: self.upd_info())
        ttk.Button(f, text="⚡ Run PCR Simulation", command=self.run).pack(pady=5)
        
        # Replace Text with SnapGeneViewer
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
            
            # Preview binding on Template sequence
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
                    # Check if swapping would work
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

            # Apply PCR specific filtering: Show only 2 selected primers
            UIHelper.render_gene_viewer(self.res, amp_seq, new_feats, 
                                        primers_to_show=[f_primer, r_primer], 
                                        lib_manager=self.lib)

            self.ln_l.config(text=f"Product Length: {len(amp_seq)} bp")
            ext_time = calculate_extension_time(len(amp_seq))
            self.ex_l.config(text=f"Extension Time: {ext_time//60}m {ext_time%60}s")

        except Exception as e:
            messagebox.showerror("Error", f"PCR Simulation failed: {str(e)}\n{traceback.format_exc()}")

    def save(self):
        if hasattr(self, 'last_amp'):
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
        ttk.Button(f, text="⚡ Run HR Simulation", command=self.run).pack(pady=5)
        
        # Replace Text with SnapGeneViewer
        self.res = SnapGeneViewer(f, lib_manager=self.lib)
        self.res.pack(fill="both", expand=True)
        
        ttk.Button(f, text="Save Result to Library", command=self.save).pack(pady=5)
        UIHelper.create_legend(f).pack(side="bottom", fill="x", pady=5)

    def run(self):
        t_sel = self.target.get().split("] ")[-1] if "]" in self.target.get() else self.target.get()
        target_obj = self.lib.templates.get(t_sel)
        amp_sel = self.amp.get()
        if amp_sel.startswith("[Digest] "):
            amp_obj = self.lib.digests.get(amp_sel[len("[Digest] "):])
        else:
            amp_obj = self.lib.amplicons.get(amp_sel) or self.lib.digests.get(amp_sel)

        if target_obj and amp_obj:
            final_seq, feats, s, e = simulate_hr_logic_with_features(target_obj, amp_obj)
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
            else:
                messagebox.showerror("Error", "No homologous recombination sites found (min 10bp match at ends).")

    def save(self):
        if hasattr(self, 'last_hr'):
            n = f"HR_{self.target.get()[:10]}_{self.amp.get()[:10]}"
            self.last_hr.name = n
            self.lib.recombinants[n] = self.last_hr
            self.lib.save()
            self.refresh_cb()
            messagebox.showinfo("Success", f"Recombinant '{n}' saved to library.")
