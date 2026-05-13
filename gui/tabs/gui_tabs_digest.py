"""제한효소 처리(Restriction Digest) 탭
- circular plasmid 또는 linear template를 source로 받음
- Biopython REBASE의 enzyme을 검색·다중 선택
- cut site 검색 → fragment 생성 → 선택 저장
- 저장된 fragment는 SequenceItem(category="digest")로 라이브러리에 들어가
  HR Simulation 등 다른 탭의 target으로 사용 가능
"""
import tkinter as tk
from tkinter import ttk, messagebox
from core.models import SequenceItem
from gui.gui_components import UIHelper, SearchableCombobox, SnapGeneViewer
from core.restriction import (
    get_all_enzyme_names, find_cut_sites, digest,
    inherit_features, enzyme_info, annotate_cut_sites_on_source,
)


class DigestTab:
    def __init__(self, parent, lib, refresh_callback):
        self.parent, self.lib, self.refresh_cb = parent, lib, refresh_callback
        self.last_fragments = []
        self.current_source_obj = None
        self.current_enzymes = []
        self.current_sites = []
        self.all_enzymes = get_all_enzyme_names()
        self.setup_ui()

    # --- UI 구성 ---
    def setup_ui(self):
        f = ttk.Frame(self.parent, padding=10); f.pack(fill="both", expand=True)

        # 1. Source selector
        sf = ttk.LabelFrame(f, text="Source (Plasmid/Template)", padding=8); sf.pack(fill="x")
        r0 = ttk.Frame(sf); r0.pack(fill="x")
        ttk.Label(r0, text="Source:", width=8).pack(side="left")
        self.src = SearchableCombobox(r0)
        self.src.pack(side="left", fill="x", expand=True, padx=5)
        self.src.bind("<<ComboboxSelected>>", lambda e: self.on_source_change())
        self.src_info = ttk.Label(sf, text="(no source selected)",
                                  font=("TkDefaultFont", 8, "italic"), foreground="gray")
        self.src_info.pack(anchor="w", pady=2)

        # 2. Enzyme selection (검색 + 좌우 listbox)
        ef = ttk.LabelFrame(f, text="Restriction Enzymes (REBASE)", padding=8); ef.pack(fill="x", pady=5)
        rs = ttk.Frame(ef); rs.pack(fill="x")
        ttk.Label(rs, text="🔍", width=2).pack(side="left")
        self.enz_search = ttk.Entry(rs)
        self.enz_search.pack(side="left", fill="x", expand=True, padx=4)
        self.enz_search.bind("<KeyRelease>", lambda e: self.refresh_enzyme_list())
        ttk.Label(rs, text=f"(total: {len(self.all_enzymes)})",
                  font=("TkDefaultFont", 8), foreground="gray").pack(side="left", padx=4)

        lr = ttk.Frame(ef); lr.pack(fill="x", pady=4)
        # available
        lf = ttk.Frame(lr); lf.pack(side="left", fill="both", expand=True)
        ttk.Label(lf, text="Available", font=("TkDefaultFont", 8, "bold")).pack(anchor="w")
        lb_fr = ttk.Frame(lf); lb_fr.pack(fill="both", expand=True)
        self.enz_list = tk.Listbox(lb_fr, height=8, selectmode="extended", exportselection=False,
                                   font=("Courier New", 9))
        sb1 = ttk.Scrollbar(lb_fr, orient="vertical", command=self.enz_list.yview)
        self.enz_list.config(yscrollcommand=sb1.set)
        self.enz_list.pack(side="left", fill="both", expand=True)
        sb1.pack(side="right", fill="y")
        self.enz_list.bind("<Double-Button-1>", lambda e: self.add_enzyme())
        self.enz_list.bind("<<ListboxSelect>>", lambda e: self.show_enzyme_info())

        # buttons
        bf = ttk.Frame(lr); bf.pack(side="left", padx=6)
        ttk.Button(bf, text="→ Add", command=self.add_enzyme, width=8).pack(pady=2)
        ttk.Button(bf, text="← Remove", command=self.remove_enzyme, width=8).pack(pady=2)
        ttk.Button(bf, text="Clear", command=lambda: self.sel_list.delete(0, tk.END), width=8).pack(pady=2)

        # selected
        rf = ttk.Frame(lr); rf.pack(side="left", fill="both", expand=True)
        ttk.Label(rf, text="Selected", font=("TkDefaultFont", 8, "bold")).pack(anchor="w")
        lb_fr2 = ttk.Frame(rf); lb_fr2.pack(fill="both", expand=True)
        self.sel_list = tk.Listbox(lb_fr2, height=8, selectmode="extended", exportselection=False,
                                   font=("Courier New", 9))
        sb2 = ttk.Scrollbar(lb_fr2, orient="vertical", command=self.sel_list.yview)
        self.sel_list.config(yscrollcommand=sb2.set)
        self.sel_list.pack(side="left", fill="both", expand=True)
        sb2.pack(side="right", fill="y")
        self.sel_list.bind("<Double-Button-1>", lambda e: self.remove_enzyme())

        self.enz_info = ttk.Label(ef, text="", font=("Courier New", 8), foreground="darkblue")
        self.enz_info.pack(anchor="w", pady=2)

        # 3. Run
        run_f = ttk.Frame(f); run_f.pack(fill="x", pady=5)
        ttk.Button(run_f, text="✂️ Run Digestion", command=self.run).pack(side="left")
        ttk.Button(run_f, text="🔍 Preview Cut Sites on Source",
                   command=self.preview_sites).pack(side="left", padx=10)
        self.summary_l = ttk.Label(run_f, text="", font=("TkDefaultFont", 9, "bold"))
        self.summary_l.pack(side="left", padx=15)

        # 4. Fragment table
        ft_lf = ttk.LabelFrame(f, text="Fragments", padding=4); ft_lf.pack(fill="x", pady=4)
        fr_t, self.frag_tree = UIHelper.create_scrolled_tree(
            ft_lf,
            ("idx", "len", "cs", "ce", "ends"),
            ("#", "Length(bp)", "Cut Start", "Cut End", "Cut by (5'→3')")
        )
        fr_t.pack(fill="x")
        self.frag_tree.column("idx", width=40)
        self.frag_tree.column("len", width=80)
        self.frag_tree.column("cs", width=80)
        self.frag_tree.column("ce", width=80)
        self.frag_tree.column("ends", width=200)
        self.frag_tree.bind("<<TreeviewSelect>>", lambda e: self.on_frag_sel())

        # 5. Detail (sequence with annotations)
        det_lf = ttk.LabelFrame(f, text="Selected Fragment Sequence", padding=4)
        det_lf.pack(fill="both", expand=True, pady=4)
        
        # Replace Text with SnapGeneViewer
        self.det = SnapGeneViewer(det_lf, lib_manager=self.lib)
        self.det.pack(fill="both", expand=True)

        # 6. Save
        sv = ttk.Frame(f); sv.pack(fill="x", pady=4)
        ttk.Label(sv, text="Save name:").pack(side="left")
        self.save_name = ttk.Entry(sv, width=40)
        self.save_name.pack(side="left", padx=4)
        ttk.Button(sv, text="💾 Save Selected", command=self.save_selected).pack(side="left", padx=4)
        ttk.Button(sv, text="💾 Save All Fragments", command=self.save_all).pack(side="left", padx=4)

        UIHelper.create_legend(f).pack(side="bottom", fill="x", pady=4)

        self.refresh_enzyme_list()

    # --- Source 콤보박스 갱신 (외부에서 호출) ---
    def refresh_sources(self):
        items = []
        for k, v in self.lib.templates.items():
            items.append(f"[{v.topology}] {k}")
        # Amplicon 추가
        for k, v in self.lib.amplicons.items():
            items.append(f"[Amplicon:{v.topology}] {k}")
        # 이미 digest 된 것도 다시 자를 수 있게
        for k, v in self.lib.digests.items():
            items.append(f"[Digest:{v.topology}] {k}")
        self.src['values'] = items

    def on_source_change(self):
        obj = self._get_source_obj()
        if obj:
            self.src_info.config(
                text=f"Length: {len(obj.sequence)} bp  |  Topology: {obj.topology}"
                     f"  |  Features: {len(obj.features)}",
                foreground="black"
            )
        else:
            self.src_info.config(text="(no source selected)", foreground="gray")

    def _get_source_obj(self):
        s = self.src.get()
        if not s: return None
        # "[Linear] name", "[Amplicon:Linear] name" 또는 "[Digest:Linear] name" 형식
        name = s.split("] ", 1)[-1] if "] " in s else s
        return self.lib.templates.get(name) or self.lib.amplicons.get(name) or self.lib.digests.get(name)

    # --- Enzyme listbox 관리 ---
    def refresh_enzyme_list(self):
        q = self.enz_search.get().strip().upper()
        self.enz_list.delete(0, tk.END)
        for e in self.all_enzymes:
            if not q or q in e.upper():
                self.enz_list.insert(tk.END, e)

    def show_enzyme_info(self):
        sel = self.enz_list.curselection()
        if not sel: return
        ename = self.enz_list.get(sel[0])
        self.enz_info.config(text=enzyme_info(ename))

    def add_enzyme(self):
        existing = set(self.sel_list.get(0, tk.END))
        for i in self.enz_list.curselection():
            ename = self.enz_list.get(i)
            if ename not in existing:
                self.sel_list.insert(tk.END, ename)
                existing.add(ename)

    def remove_enzyme(self):
        for i in reversed(self.sel_list.curselection()):
            self.sel_list.delete(i)

    # --- Cut site preview ---
    def preview_sites(self):
        src = self._get_source_obj()
        if not src:
            messagebox.showwarning("Warning", "Source를 먼저 선택하세요.")
            return
        enzymes = list(self.sel_list.get(0, tk.END))
        if not enzymes:
            messagebox.showwarning("Warning", "효소를 1개 이상 선택하세요.")
            return
        sites = find_cut_sites(src.sequence, src.topology, enzymes)
        if not sites:
            messagebox.showinfo("No cuts", "선택한 효소들이 이 서열을 자르지 않습니다.")
            return
        # 원본 서열에 cut site를 feature로 추가하여 표시
        feats = list(src.features) + annotate_cut_sites_on_source(src.sequence, sites)
        UIHelper.render_gene_viewer(self.det, src.sequence, feats, lib_manager=self.lib)
        summary = f"{len(sites)} cut(s): " + ", ".join(
            f"{s['enzyme']}@{s['pos']}" for s in sites
        )
        self.summary_l.config(text=summary, foreground="darkgreen")

    # --- Run digestion ---
    def run(self):
        src = self._get_source_obj()
        if not src:
            messagebox.showwarning("Warning", "Source를 먼저 선택하세요.")
            return
        enzymes = list(self.sel_list.get(0, tk.END))
        if not enzymes:
            messagebox.showwarning("Warning", "효소를 1개 이상 선택하세요.")
            return

        sites = find_cut_sites(src.sequence, src.topology, enzymes)
        if not sites:
            messagebox.showinfo("No cuts", "선택한 효소들이 이 서열을 자르지 않습니다.")
            self.summary_l.config(text="No cut sites found.", foreground="red")
            return

        fragments = digest(src.sequence, src.topology, sites)
        for frag in fragments:
            frag["features"] = inherit_features(
                src.features, frag, src.topology, len(src.sequence)
            )

        self.last_fragments = fragments
        self.current_source_obj = src
        self.current_enzymes = enzymes
        self.current_sites = sites

        # fragment 테이블 갱신
        self.frag_tree.delete(*self.frag_tree.get_children())
        for i, f in enumerate(fragments):
            ends = f"{f.get('start_enzyme','?')}  ↔  {f.get('end_enzyme','?')}"
            self.frag_tree.insert("", "end", iid=str(i), values=(
                i + 1, len(f["sequence"]), f["cut_start"], f["cut_end"], ends
            ))

        self.summary_l.config(
            text=f"{len(sites)} cut(s) → {len(fragments)} fragment(s).",
            foreground="darkgreen"
        )
        # 첫 fragment 자동 선택
        if fragments:
            self.frag_tree.selection_set("0")
            self.on_frag_sel()

    def on_frag_sel(self):
        sel = self.frag_tree.selection()
        if not sel: return
        i = int(sel[0])
        if i >= len(self.last_fragments): return
        frag = self.last_fragments[i]
        UIHelper.render_gene_viewer(self.det, frag["sequence"], frag["features"], lib_manager=self.lib)
        # 자동 이름 제안
        self.save_name.delete(0, tk.END)
        self.save_name.insert(0, self._make_name(i))

    def _make_name(self, frag_idx):
        if not self.current_source_obj: return ""
        enz_str = "+".join(self.current_enzymes) if len(self.current_enzymes) > 1 else self.current_enzymes[0]
        base = f"{enz_str} cut {self.current_source_obj.name}"
        if len(self.last_fragments) > 1:
            base += f" #{frag_idx + 1}"
        return base

    # --- 저장 ---
    def _save_one(self, idx, override_name=None):
        if idx >= len(self.last_fragments): return False
        frag = self.last_fragments[idx]
        name = (override_name or self._make_name(idx)).strip()
        if not name: return False
        # 이름 중복 회피
        original = name
        ctr = 2
        while name in self.lib.digests:
            name = f"{original} ({ctr})"
            ctr += 1
        item = SequenceItem(
            name=name,
            sequence=frag["sequence"],
            category="digest",
            topology="Linear",
            features=frag["features"],
            template_name=self.current_source_obj.name,
            enzymes=",".join(self.current_enzymes),
        )
        self.lib.digests[name] = item
        return True

    def save_selected(self):
        sel = self.frag_tree.selection()
        if not sel:
            messagebox.showwarning("Warning", "저장할 fragment를 선택하세요.")
            return
        # 1개만 선택 시 사용자 입력 이름 사용, 여러 개면 자동 이름
        if len(sel) == 1:
            name = self.save_name.get().strip()
            if self._save_one(int(sel[0]), override_name=name):
                self.lib.save(); self.refresh_cb()
                messagebox.showinfo("Saved", f"'{name}' 저장 완료.\nLibrary > Digest Product 탭에서 확인.")
        else:
            ok = 0
            for sid in sel:
                if self._save_one(int(sid)): ok += 1
            self.lib.save(); self.refresh_cb()
            messagebox.showinfo("Saved", f"{ok}개 fragment 저장 완료.")

    def save_all(self):
        if not self.last_fragments:
            messagebox.showwarning("Warning", "먼저 Run Digestion을 실행하세요.")
            return
        ok = 0
        for i in range(len(self.last_fragments)):
            if self._save_one(i): ok += 1
        self.lib.save(); self.refresh_cb()
        messagebox.showinfo("Saved", f"{ok}개 fragment 저장 완료.")
