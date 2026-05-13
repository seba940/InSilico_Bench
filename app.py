import tkinter as tk
from tkinter import ttk
from core.manager import LibraryManager
from gui.tabs.gui_tabs_lib import PrimerTab, TemplateTab, MarkerTab, TagTab
from gui.tabs.gui_tabs_tools import SpeciesTab, AnnotationRefTab, SeqViewTab, PrimerValidationTab
from gui.tabs.gui_tabs_sim import PCRTab, HRTab, BatchPCRTab
from gui.tabs.gui_tabs_digest import DigestTab
from gui.tabs.gui_tabs_ligation import LigationTab
from gui.tabs.gui_tabs_genelookup import GeneLookupTab
from gui.tabs.gui_tabs_hrdesign import HRDesignTab
from gui.gui_components import ScrollableFrame


class YeastApp:
    def __init__(self, root):
        self.root = root
        self.root.title("InSilico_Bench - Pro v14")
        self.root.geometry("1400x850")
        self.lib = LibraryManager()

        self.main_container = ScrollableFrame(root)
        self.main_container.pack(expand=True, fill="both")
        content_f = self.main_container.scrollable_frame

        self.main_nb = ttk.Notebook(content_f)
        self.main_nb.pack(expand=True, fill="both", padx=10, pady=10)

        # --- Library Management tab ---
        self.tab_lib = ttk.Frame(self.main_nb)
        self.main_nb.add(self.tab_lib, text=" \U0001f4c1 Library Management ")

        self.lib_nb = ttk.Notebook(self.tab_lib)
        self.lib_nb.pack(expand=True, fill="both", padx=5, pady=5)

        self.f_p   = ttk.Frame(self.lib_nb)
        self.f_t   = ttk.Frame(self.lib_nb)
        self.f_a   = ttk.Frame(self.lib_nb)
        self.f_r   = ttk.Frame(self.lib_nb)
        self.f_d   = ttk.Frame(self.lib_nb)
        self.f_lig = ttk.Frame(self.lib_nb)
        self.f_m   = ttk.Frame(self.lib_nb)
        self.f_tag = ttk.Frame(self.lib_nb)
        self.f_ann = ttk.Frame(self.lib_nb)
        self.f_sp  = ttk.Frame(self.lib_nb)

        self.lib_nb.add(self.f_p,   text=" Primer ")
        self.lib_nb.add(self.f_t,   text=" Template ")
        self.lib_nb.add(self.f_a,   text=" Amplicon ")
        self.lib_nb.add(self.f_r,   text=" HR Results ")
        self.lib_nb.add(self.f_d,   text=" Digest Product ")
        self.lib_nb.add(self.f_lig, text=" Ligation Product ")
        self.lib_nb.add(self.f_m,   text=" Selection Marker ")
        self.lib_nb.add(self.f_tag, text=" Epitope Tag ")
        self.lib_nb.add(self.f_ann, text=" Annotation Reference ")
        self.lib_nb.add(self.f_sp,  text=" Species Manager ")

        self.p_tab       = PrimerTab(self.f_p,   self.lib, self.refresh_all)
        self.t_tab       = TemplateTab(self.f_t, self.lib, self.refresh_all)
        self.a_tab       = SeqViewTab(self.f_a,   self.lib, self.refresh_all, self.lib.amplicons,   "amplicon")
        self.r_tab       = SeqViewTab(self.f_r,   self.lib, self.refresh_all, self.lib.recombinants, "recombinant")
        self.d_tab       = SeqViewTab(self.f_d,   self.lib, self.refresh_all, self.lib.digests,     "digest")
        self.lig_lib_tab = SeqViewTab(self.f_lig, self.lib, self.refresh_all, self.lib.ligations,   "ligation")
        self.m_tab       = MarkerTab(self.f_m,   self.lib, self.refresh_all)
        self.tag_tab     = TagTab(self.f_tag,    self.lib, self.refresh_all)
        self.ann_tab     = AnnotationRefTab(self.f_ann, self.lib, self.refresh_all)
        self.sp_tab      = SpeciesTab(self.f_sp, self.lib, self.refresh_all)

        # --- Simulation / tool tabs ---
        self.f_pcr     = ttk.Frame(self.main_nb)
        self.f_bpcr    = ttk.Frame(self.main_nb)
        self.f_hr      = ttk.Frame(self.main_nb)
        self.f_dig     = ttk.Frame(self.main_nb)
        self.f_lig_tab = ttk.Frame(self.main_nb)
        self.f_val     = ttk.Frame(self.main_nb)
        self.f_gl      = ttk.Frame(self.main_nb)
        self.f_hrd     = ttk.Frame(self.main_nb)

        self.main_nb.add(self.f_pcr,     text=" \U0001f9ec PCR Execution ")
        self.main_nb.add(self.f_bpcr,    text=" \U0001f4cb Batch PCR ")
        self.main_nb.add(self.f_hr,      text=" \U0001f504 HR Simulation ")
        self.main_nb.add(self.f_hrd,     text=" \U0001f9ec HR Primer Design ")
        self.main_nb.add(self.f_dig,     text=" ✂️ Restriction Digest ")
        self.main_nb.add(self.f_lig_tab, text=" \U0001f517 Ligation ")
        self.main_nb.add(self.f_val,     text=" ✅ Primer Validation ")
        self.main_nb.add(self.f_gl,      text=" \U0001f9ec Gene Lookup ")

        self.pcr_tab  = PCRTab(self.f_pcr,  self.lib, self.refresh_all)
        self.bpcr_tab = BatchPCRTab(self.f_bpcr, self.lib, self.refresh_all)
        self.hr_tab   = HRTab(self.f_hr,    self.lib, self.refresh_all)
        self.hrd_tab  = HRDesignTab(self.f_hrd, self.lib, self.refresh_all)
        self.dig_tab  = DigestTab(self.f_dig, self.lib, self.refresh_all)
        self.lig_tab  = LigationTab(self.f_lig_tab, self.lib, self.refresh_all)
        self.val_tab  = PrimerValidationTab(self.f_val, self.lib, self.refresh_all)
        self.gl_tab   = GeneLookupTab(self.f_gl,  self.lib, self.refresh_all)

        self.main_nb.bind("<<NotebookTabChanged>>", lambda e: self.refresh_all())
        self.root.bind("<Control-f>", self.show_find_dialog)
        self.root.bind("<Control-F>", self.show_find_dialog)
        self.root.bind("<Control-z>", self.undo)
        self.root.bind("<Control-Z>", self.undo)
        self.refresh_all()

    def show_find_dialog(self, event=None):
        curr = self.main_nb.index("current")
        txt, f_ed = None, None
        if curr == 0:
            sub = self.lib_nb.index("current")
            if sub == 1:  txt, f_ed = self.t_tab.ts,       self.t_tab.feat_ed
            elif sub == 2: txt, f_ed = self.a_tab.det,      self.a_tab.feat_ed
            elif sub == 3: txt, f_ed = self.r_tab.det,      self.r_tab.feat_ed
            elif sub == 4: txt, f_ed = self.d_tab.det,      self.d_tab.feat_ed
            elif sub == 6: txt, f_ed = self.m_tab.ms,       self.m_tab.feat_ed
            elif sub == 7: txt, f_ed = self.tag_tab.ts,     self.tag_tab.feat_ed
        elif curr == 1: txt = self.pcr_tab.res
        elif curr == 2: txt = self.hr_tab.res
        elif curr == 3: txt = self.dig_tab.det
        elif curr == 5: txt = self.val_tab.res_text
        elif curr == 6: txt = self.gl_tab.det
        if txt:
            from gui.gui_components import FindDialog
            FindDialog(self.root, txt, feature_editor=f_ed)

    def undo(self, event=None):
        if self.lib.undo():
            self.refresh_all()
            from tkinter import messagebox
            messagebox.showinfo("Undo", "마지막 라이브러리 변경이 취소되었습니다.")
        else:
            from tkinter import messagebox
            messagebox.showinfo("Undo", "더 이상 되돌릴 수 없습니다.")

    def refresh_all(self):
        self.p_tab.filter_p(self.p_tab.search.get())
        self.t_tab.filter(self.t_tab.search.get())
        self.a_tab.filter(self.a_tab.search.get())
        self.r_tab.filter(self.r_tab.search.get())
        self.d_tab.filter(self.d_tab.search.get())
        self.lig_lib_tab.filter(self.lig_lib_tab.search.get())
        self.m_tab.filter(self.m_tab.search.get())
        self.tag_tab.filter(self.tag_tab.search.get())
        self.ann_tab.filter(self.ann_tab.search.get())
        self.sp_tab.filter(self.sp_tab.search.get())
        self.val_tab.refresh_lists()

        sl = sorted(list(self.lib.species.keys()))
        self.t_tab.sp_c['values'] = sl

        ml = [""] + sorted(list(self.lib.markers.keys()))
        tl = [""] + sorted(list(self.lib.tags.keys()))
        for t in [self.a_tab, self.r_tab, self.d_tab, self.lig_lib_tab]:
            t.ae_m['values'] = ml
            t.ae_t['values'] = tl

        # PCR template: template + digest + HR result + ligation product
        tl_all = [f"[{v.topology}] {k}"          for k, v in self.lib.templates.items()]
        dl_all = [f"[Digest:{v.topology}] {k}"   for k, v in self.lib.digests.items()]
        rl_all = [f"[HR:{v.topology}] {k}"       for k, v in self.lib.recombinants.items()]
        ll_all = [f"[Ligation:{v.topology}] {k}" for k, v in self.lib.ligations.items()]
        self.pcr_tab.tmp['values'] = tl_all + dl_all + rl_all + ll_all

        # HR target: Genome-kind templates only
        def _kind(v):
            return v.kind or ("Genome" if v.topology == "Linear" else "Plasmid")
        hr_target_list = [f"[{v.topology}] {k}"
                          for k, v in self.lib.templates.items() if _kind(v) == "Genome"]
        self.hr_tab.target['values'] = hr_target_list

        # HR insert: amplicon + digest + ligation
        amp_list = sorted(list(self.lib.amplicons.keys()))
        dig_list = [f"[Digest] {k}"   for k in sorted(self.lib.digests.keys())]
        lig_list = [f"[Ligation] {k}" for k in sorted(self.lib.ligations.keys())]
        self.hr_tab.amp['values'] = amp_list + dig_list + lig_list

        self.pcr_tab.fwd['values'] = sorted(list(self.lib.fwd_primers.keys()))
        self.pcr_tab.rev['values'] = sorted(list(self.lib.rev_primers.keys()))
        self.bpcr_tab.refresh_combos()

        # HR colony PCR primer lists
        hr_fwd = sorted(list(self.lib.fwd_primers.keys()))
        hr_rev = sorted(list(self.lib.rev_primers.keys()))
        self.hr_tab.col_fwd['values'] = hr_fwd
        self.hr_tab.col_rev['values'] = hr_rev

        # source combos
        self.dig_tab.refresh_sources()
        self.lig_tab.refresh_sources()
        self.gl_tab.refresh_sources()
        self.hrd_tab.refresh_sources()
