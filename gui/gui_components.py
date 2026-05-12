import tkinter as tk
from tkinter import ttk, messagebox
try:
    from Bio.SeqUtils import MeltingTemp as mt
    from Bio.SeqUtils import gc_fraction
except ImportError:
    mt = None
    gc_fraction = None

class Tooltip:
    """실시간으로 내용이 변하는 지능형 말풍선"""
    def __init__(self, widget):
        self.widget = widget
        self.tip_window = None
        self.label = None

    def show(self, text, x, y):
        if self.tip_window:
            # 이미 창이 있다면 텍스트와 위치만 갱신
            self.label.config(text=text)
            self.tip_window.wm_geometry(f"+{x+15}+{y+10}")
            return
        
        self.tip_window = tw = tk.Toplevel(self.widget)
        tw.wm_overrideredirect(True)
        tw.wm_geometry(f"+{x+15}+{y+10}")
        self.label = tk.Label(tw, text=text, justify=tk.LEFT, background="#ffffe0", 
                              relief=tk.SOLID, borderwidth=1, font=("TkDefaultFont", "8", "normal"))
        self.label.pack()

    def hide(self):
        tw = self.tip_window
        self.tip_window = None
        if tw: tw.destroy()

class ScrollableFrame(ttk.Frame):
    def __init__(self, container, *args, **kwargs):
        super().__init__(container, *args, **kwargs)
        self.canvas = tk.Canvas(self, highlightthickness=0)
        self.scrollbar_v = ttk.Scrollbar(self, orient="vertical", command=self.canvas.yview)
        self.scrollbar_h = ttk.Scrollbar(self, orient="horizontal", command=self.canvas.xview)
        self.scrollable_frame = ttk.Frame(self.canvas)
        self.scrollable_frame.bind("<Configure>", lambda e: self.canvas.configure(scrollregion=self.canvas.bbox("all")))
        self.canvas_window = self.canvas.create_window((0, 0), window=self.scrollable_frame, anchor="nw")
        self.canvas.bind('<Configure>', lambda e: self.canvas.itemconfig(self.canvas_window, width=e.width))
        self.canvas.configure(yscrollcommand=self.scrollbar_v.set, xscrollcommand=self.scrollbar_h.set)
        self.canvas.bind_all("<MouseWheel>", lambda e: self.canvas.yview_scroll(int(-1*(e.delta/120)), "units"))
        self.canvas.pack(side="left", fill="both", expand=True)
        self.scrollbar_v.pack(side="right", fill="y")
        self.scrollbar_h.pack(side="bottom", fill="x")

class UIHelper:
    COLOR_MAP = {
        "Misc": "#F5DEB3", "Homology Arm": "#E0E0E0", "Marker": "#FFB2B2",
        "Tag": "#E6B2FF", "Terminator": "#FFFFB2", "Promoter": "#B2FFB2",
        "CDS": "#B2D8FF", "Primer": "#FFA040"
    }

    @staticmethod
    def create_scrolled_tree(parent, columns, headings, lib=None, table_id=None):
        frame = ttk.Frame(parent)
        tree = ttk.Treeview(frame, columns=columns, show="headings", selectmode="extended", height=5)
        vsb = ttk.Scrollbar(frame, orient="vertical", command=tree.yview)
        tree.configure(yscrollcommand=vsb.set)
        
        # 최적화된 기본 넓이 설정 (데이터 종류에 따라 다르게 설정 가능)
        default_widths = {
            "N": 150, "Name": 150, "Usage": 80, "U": 80, 
            "Tm": 60, "Tm_B": 70, "Tm_T": 70, "S": 250, "Seq": 250, "L": 60, "Len": 60,
            "Topology": 80, "T": 80, "Kind": 80, "K": 80, "Species": 120
        }

        for col, head in zip(columns, headings):
            tree.heading(col, text=head)
            # 저장된 넓이가 있으면 로드, 없으면 기본값
            w = 100
            if lib and table_id and "col_widths" in lib.settings:
                w = lib.settings["col_widths"].get(f"{table_id}_{col}", default_widths.get(col, 100))
            else:
                w = default_widths.get(col, 100)
            tree.column(col, width=w, minwidth=50, stretch=True)

        if lib and table_id:
            def save_widths(event):
                if tree.identify_region(event.x, event.y) == "separator":
                    if "col_widths" not in lib.settings: lib.settings["col_widths"] = {}
                    for c in columns:
                        lib.settings["col_widths"][f"{table_id}_{c}"] = tree.column(c, 'width')
                    lib.save()
            tree.bind("<ButtonRelease-1>", save_widths)

        tree.pack(side="left", fill="both", expand=True); vsb.pack(side="right", fill="y")
        return frame, tree

    @staticmethod
    def add_search_bar(parent, filter_cmd, dup_cmd):
        f = ttk.Frame(parent); f.pack(fill="x", pady=1)
        ttk.Label(f, text="🔍 Search:", font=("TkDefaultFont", 8)).pack(side="left", padx=2)
        ent = ttk.Entry(f, font=("TkDefaultFont", 8)); ent.pack(side="left", fill="x", expand=True, padx=2)
        ent.bind("<KeyRelease>", lambda e: filter_cmd(ent.get()))
        if dup_cmd: ttk.Button(f, text="✨ Deduplicate", command=dup_cmd, width=10).pack(side="right", padx=2)
        return ent

    @staticmethod
    def setup_inline_edit(tree, columns, values_map, on_save):
        sel = tree.selection()
        if not sel: return []
        item_id = sel[0]; widgets = []
        for col_id, choices in columns.items():
            rect = tree.bbox(item_id, col_id)
            if not rect: continue
            w = ttk.Combobox(tree, values=choices, state="readonly", font=("TkDefaultFont", 8)) if choices else ttk.Entry(tree, font=("TkDefaultFont", 8))
            if not choices: w.insert(0, tree.item(item_id)['values'][values_map[col_id]])
            else: w.set(tree.item(item_id)['values'][values_map[col_id]])
            w.place(x=rect[0], y=rect[1], width=rect[2], height=rect[3])
            w.bind("<Return>", lambda e: on_save(item_id))
            widgets.append(w)
        return widgets

    @classmethod
    def render_gene_viewer(cls, viewer, sequence, features=None, primers_to_show=None, lib_manager=None):
        """
        Unified helper to render sequence in SnapGeneViewer.
        - viewer: SnapGeneViewer instance
        - sequence: string
        - features: list of feature dicts (optional)
        - primers_to_show: 
            - None: Auto-detect all primers from global library (default)
            - list of primer objects: show only these specific primers
            - []: show no primers
        - lib_manager: PrimerLibrary manager for global features and primers
        """
        if not isinstance(viewer, SnapGeneViewer):
            # Fallback for legacy tk.Text if still used somewhere
            cls.render_annotations(viewer, sequence, features, lib_manager, primers_to_show)
            return

        cls.render_annotations(viewer, sequence, features, lib_manager, primers_to_show)

    @classmethod
    def render_annotations(cls, text_widget, sequence, features, lib_manager=None, primers_to_show=None):
        # ── SnapGeneViewer fast-path ─────────────────────────────────
        all_features = list(features) if features else []
        
        # 1. Global features (auto-annotation)
        if lib_manager and hasattr(lib_manager, 'global_features'):
            from core.utils import auto_annotate
            auto_feats = auto_annotate(sequence, lib_manager.global_features)
            seen_pos = set((f['start'], f['end'], f.get('label')) for f in all_features)
            for af in auto_feats:
                if (af['start'], af['end'], af.get('label')) not in seen_pos:
                    all_features.append(af)

        # 2. Primers (binding site detection)
        from core.utils import find_primer_bindings
        if primers_to_show is not None:
            # Show only specific primers (e.g. PCR execution tab)
            fwd_list = [p for p in primers_to_show if getattr(p, 'p_type', '') == 'Fwd']
            rev_list = [p for p in primers_to_show if getattr(p, 'p_type', '') == 'Rev']
            
            # If they are just names, resolve them
            if fwd_list and isinstance(fwd_list[0], str) and lib_manager:
                fwd_list = [lib_manager.fwd_primers.get(p) for p in fwd_list if lib_manager.fwd_primers.get(p)]
            if rev_list and isinstance(rev_list[0], str) and lib_manager:
                rev_list = [lib_manager.rev_primers.get(p) for p in rev_list if lib_manager.rev_primers.get(p)]

            primer_feats = find_primer_bindings(sequence, fwd_list, rev_list)
            all_features.extend(primer_feats)
        elif lib_manager and hasattr(lib_manager, 'fwd_primers'):
            # Default: Show all primers from global library
            primer_feats = find_primer_bindings(
                sequence,
                lib_manager.fwd_primers,
                lib_manager.rev_primers)
            seen_primer = set(
                (f['start'], f['end'], f.get('strand', 1), f.get('label'))
                for f in all_features if f.get('type') == 'Primer')
            for pf in primer_feats:
                key = (pf['start'], pf['end'], pf['strand'], pf.get('label'))
                if key not in seen_primer:
                    seen_primer.add(key)
                    all_features.append(pf)

        if isinstance(text_widget, SnapGeneViewer):
            if text_widget.lib is None:
                text_widget.lib = lib_manager
                # Re-load width if it was just attached
                if lib_manager and "viewer_row_width" in lib_manager.settings:
                    text_widget._row_width = lib_manager.settings["viewer_row_width"]
                    text_widget._rw_var.set(text_widget._row_width)
            text_widget.render(sequence, all_features)
            return
        
        # ── legacy tk.Text path ──────────────────────────────────────
        text_widget.config(state="normal", font=("Courier New", 9))
        text_widget.delete("1.0", tk.END); text_widget.insert(tk.END, sequence)
        
        text_widget.feature_data = all_features
        if not hasattr(text_widget, 'tooltip'):
            text_widget.tooltip = Tooltip(text_widget)
            def on_motion(event):
                try:
                    index = text_widget.index(f"@{event.x},{event.y}")
                    char_idx = int(index.split('.')[1])
                    found = None
                    # 짧은 특징 우선 탐색
                    for f in sorted(text_widget.feature_data, key=lambda x: (x['end']-x['start'])):
                        if f['start'] <= char_idx < f['end']:
                            found = f; break
                    if found:
                        text_widget.tooltip.show(f"[{found.get('type','Misc')}] {found['label']}", event.x_root, event.y_root)
                    else: text_widget.tooltip.hide()
                except: pass
            text_widget.bind("<Motion>", on_motion)
            text_widget.bind("<Leave>", lambda e: text_widget.tooltip.hide())

        for f_type, color in cls.COLOR_MAP.items():
            text_widget.tag_configure(f_type, background=color, foreground="black")
            text_widget.tag_raise(f_type)
        text_widget.tag_configure("search_hit", background="orange", foreground="black")
        text_widget.tag_raise("search_hit")

        for f in sorted(all_features, key=lambda x: (x['end'] - x['start']), reverse=True):
            tag = f.get('type', 'Misc')
            if tag not in cls.COLOR_MAP: tag = "Misc"
            try: text_widget.tag_add(tag, f"1.{f['start']}", f"1.{f['end']}")
            except: pass

    @classmethod
    def create_legend(cls, parent_frame, orient="horizontal"):
        l_frame = ttk.LabelFrame(parent_frame, text="🎨 Legend", padding=1)
        if orient == "horizontal":
            row_f = ttk.Frame(l_frame); row_f.pack(fill="x")
            for i, (f_type, color) in enumerate(cls.COLOR_MAP.items()):
                item_f = ttk.Frame(row_f); item_f.pack(side="left", padx=3)
                tk.Label(item_f, bg=color, width=1, height=0, relief="ridge").pack(side="left")
                ttk.Label(item_f, text=f_type, font=("TkDefaultFont", 7, "bold")).pack(side="left", padx=1)
        else: # vertical
            for f_type, color in cls.COLOR_MAP.items():
                item_f = ttk.Frame(l_frame); item_f.pack(fill="x", padx=3, pady=1)
                tk.Label(item_f, bg=color, width=1, height=0, relief="ridge").pack(side="left")
                ttk.Label(item_f, text=f_type, font=("TkDefaultFont", 7, "bold")).pack(side="left", padx=1)
        return l_frame

class SearchableCombobox(ttk.Frame):
    """Combobox with live substring-search filter.

    Drop-in replacement for ``ttk.Combobox(state="readonly")``.
    Supports:
      • widget['values'] = list  / widget.configure(values=list)
      • widget.get() / widget.set(value)
      • widget.bind("<<ComboboxSelected>>", callback)
      • Typing a partial name filters the dropdown in real-time
      • Clicking ▼ shows the full list
      • ↓ key moves focus into the list; Esc closes it
    """

    def __init__(self, parent, **kwargs):
        self._values = list(kwargs.pop("values", []))
        kwargs.pop("state", None)   # ignore readonly/normal – we manage it
        super().__init__(parent, **kwargs)

        self._var = tk.StringVar()
        self._dropdown_win = None
        self._listbox = None

        self._entry = ttk.Entry(self, textvariable=self._var)
        self._entry.pack(side="left", fill="x", expand=True)
        self._btn = ttk.Button(self, text="▼", width=2,
                               command=self._toggle_dropdown)
        self._btn.pack(side="right")

        self._entry.bind("<KeyRelease>", self._on_key)
        self._entry.bind("<Down>",       self._focus_dropdown)
        self._entry.bind("<Escape>",     lambda e: self._hide_dropdown())
        self._entry.bind("<FocusOut>",   self._schedule_hide)

    # ── compatibility interface ───────────────────────────────────────────────

    def configure(self, cnf=None, **kw):
        if cnf:
            kw.update(cnf)
        if "values" in kw:
            self._values = list(kw.pop("values"))
        if kw:
            super().configure(**kw)

    config = configure

    def __setitem__(self, key, value):
        if key == "values":
            self._values = list(value)
        else:
            super().__setitem__(key, value)

    def get(self):
        return self._var.get()

    def set(self, value):
        self._var.set(value)

    def bind(self, sequence, func=None, add=None):
        # <<ComboboxSelected>> is generated on self (the Frame); forward the rest to Entry
        if sequence == "<<ComboboxSelected>>":
            return super().bind(sequence, func, add)
        return self._entry.bind(sequence, func, add)

    # ── internal helpers ──────────────────────────────────────────────────────

    def _get_matches(self, query=""):
        if not query:
            return list(self._values)
        q = query.lower()
        return [v for v in self._values if q in v.lower()]

    def _on_key(self, event):
        if event.keysym == "Escape":
            self._hide_dropdown(); return
        if event.keysym in ("Return", "Tab"):
            return
        matches = self._get_matches(self._var.get())
        if matches:
            self._show_dropdown(matches)
        else:
            self._hide_dropdown()

    def _toggle_dropdown(self):
        if self._dropdown_win:
            self._hide_dropdown()
        else:
            self._show_dropdown(self._values)
        self._entry.focus_set()

    def _show_dropdown(self, items):
        self._hide_dropdown()
        if not items:
            return
        self.update_idletasks()
        x = self.winfo_rootx()
        y = self.winfo_rooty() + self.winfo_height()
        w = max(self.winfo_width(), 120)

        win = tk.Toplevel(self)
        win.wm_overrideredirect(True)
        win.attributes("-topmost", True)
        win.wm_geometry(f"{w}x4+{x}+{y}")
        self._dropdown_win = win

        outer = ttk.Frame(win, relief="solid", borderwidth=1)
        outer.pack(fill="both", expand=True)

        sb = ttk.Scrollbar(outer, orient="vertical")
        h = min(10, len(items))
        lb = tk.Listbox(outer, yscrollcommand=sb.set, height=h,
                        selectmode="single", exportselection=False,
                        font=("TkDefaultFont", 9))
        sb.config(command=lb.yview)
        sb.pack(side="right", fill="y")
        lb.pack(side="left", fill="both", expand=True)
        self._listbox = lb

        for item in items:
            lb.insert(tk.END, item)

        lb.bind("<<ListboxSelect>>", self._on_lb_select)
        lb.bind("<Return>",          self._on_lb_enter)
        lb.bind("<Escape>",          lambda e: self._hide_dropdown())
        lb.bind("<FocusOut>",        self._schedule_hide)

        win.update_idletasks()
        actual_h = lb.winfo_reqheight() + 4
        win.wm_geometry(f"{w}x{actual_h}+{x}+{y}")

    def _on_lb_select(self, event=None):
        if self._listbox:
            sel = self._listbox.curselection()
            if sel:
                self._var.set(self._listbox.get(sel[0]))
                self._hide_dropdown()
                self.event_generate("<<ComboboxSelected>>")

    def _on_lb_enter(self, event=None):
        self._on_lb_select()

    def _focus_dropdown(self, event=None):
        if not self._dropdown_win:
            self._show_dropdown(self._get_matches(self._var.get()) or self._values)
        if self._listbox and self._listbox.size() > 0:
            self._listbox.focus_set()
            self._listbox.select_set(0)
            self._listbox.see(0)

    def _schedule_hide(self, event=None):
        self.after(150, self._maybe_hide)

    def _maybe_hide(self):
        try:
            focused = self.focus_get()
            if focused is self._entry:
                return
            if self._listbox and focused is self._listbox:
                return
        except Exception:
            pass
        self._hide_dropdown()

    def _hide_dropdown(self):
        win = self._dropdown_win
        self._dropdown_win = None
        self._listbox = None
        if win:
            try:
                win.destroy()
            except Exception:
                pass


class FeatureEditor(ttk.LabelFrame):
    def __init__(self, parent, title="Feature Annotation", text_widget=None):
        super().__init__(parent, text=title, padding=5)
        self.features, self.target_text = [], text_widget
        self.setup_ui()

    def setup_ui(self):
        r0 = ttk.Frame(self); r0.pack(fill="x")
        ttk.Label(r0, text="L:", font=("TkDefaultFont", 8)).pack(side="left")
        self.en = ttk.Entry(r0, width=8, font=("TkDefaultFont", 8)); self.en.pack(side="left", padx=1)
        ttk.Button(r0, text="📍Set", command=self.set_from_sel, width=4).pack(side="left", padx=1)
        ttk.Label(r0, text="S:", font=("TkDefaultFont", 8)).pack(side="left")
        self.es = ttk.Entry(r0, width=4, font=("TkDefaultFont", 8)); self.es.pack(side="left")
        ttk.Label(r0, text="E:", font=("TkDefaultFont", 8)).pack(side="left")
        self.ee = ttk.Entry(r0, width=4, font=("TkDefaultFont", 8)); self.ee.pack(side="left")
        self.et = ttk.Combobox(r0, values=list(UIHelper.COLOR_MAP.keys()), width=8, state="readonly", font=("TkDefaultFont", 8))
        self.et.set("CDS"); self.et.pack(side="left", padx=1)
        ttk.Button(r0, text="+", command=self.add_feat, width=2).pack(side="left")
        f_tree, self.tree = UIHelper.create_scrolled_tree(self, ("L", "S", "E", "T"), ("L", "S", "E", "T"))
        f_tree.pack(fill="x", pady=2); self.tree.column("L", width=80)
        btn_row = ttk.Frame(self); btn_row.pack(fill="x")
        ttk.Button(btn_row, text="Del", command=self.del_feat, width=4).pack(side="left")
        ttk.Button(btn_row, text="Clear", command=self.clear_feats, width=5).pack(side="right")

    def set_pos(self, s, e):
        self.es.delete(0, tk.END); self.es.insert(0, str(s))
        self.ee.delete(0, tk.END); self.ee.insert(0, str(e))

    def set_from_sel(self):
        if not self.target_text: return
        try:
            s_idx, e_idx = self.target_text.index("sel.first"), self.target_text.index("sel.last")
            self.set_pos(int(s_idx.split('.')[1]), int(e_idx.split('.')[1]))
        except: pass

    def add_feat(self):
        try:
            f = {"label": self.en.get().strip() or self.et.get(), "start": int(self.es.get()), "end": int(self.ee.get()), "type": self.et.get()}
            self.features.append(f); self.update_tree()
        except: pass

    def del_feat(self):
        sel = self.tree.selection()
        if sel: idx = self.tree.index(sel[0]); del self.features[idx]; self.update_tree()

    def clear_feats(self): self.features = []; self.update_tree()
    def update_tree(self):
        self.tree.delete(*self.tree.get_children())
        for f in self.features: self.tree.insert("", "end", values=(f['label'], f['start'], f['end'], f['type']))
    def get_features(self): return self.features
    def set_features(self, feats): self.features = list(feats) if feats else []; self.update_tree()

class SnapGeneViewer(ttk.Frame):
    """SnapGene-style multi-row sequence viewer with feature tracks.

    Usage:
        viewer = SnapGeneViewer(parent)
        viewer.render(sequence, features)
        viewer.show_message("some text")   # for errors / status
    """

    ROW_WIDTHS = [40, 60, 80, 100, 120]

    def __init__(self, parent, row_width=60, lib_manager=None, **kwargs):
        super().__init__(parent, **kwargs)
        import tkinter.font as tkfont
        self.lib = lib_manager
        
        # Load saved row width if available
        if self.lib and "viewer_row_width" in self.lib.settings:
            row_width = self.lib.settings["viewer_row_width"]
            
        self._row_width = row_width
        self.sequence = ""
        self.features = []
        self._hit_map = []          # list of (x1,y1,x2,y2, feat) for tooltips
        self._row_map = []          # list of (row_s, row_e, y_top, y_bot) for hit-test
        self._sel_start = -1        # drag-selection anchor (abs sequence pos)
        self._sel_end   = -1        # drag-selection cursor (abs sequence pos)
        self._dragging  = False

        # ── fonts ──────────────────────────────────────────────────────
        self._font_seq   = tkfont.Font(family="Courier New", size=9)
        self._font_label = tkfont.Font(family="Courier New", size=8)
        self._font_idx   = tkfont.Font(family="Courier New", size=8, slant="italic")
        self.CW  = self._font_seq.measure("A")       # char width
        self.CH  = self._font_seq.metrics("linespace") + 2  # seq row height
        self.FH  = self.CH - 3                       # feature box height
        self.IDX = self._font_idx.measure("000000") + 10   # index column width

        # ── toolbar ────────────────────────────────────────────────────
        tb = ttk.Frame(self)
        tb.pack(fill="x")
        ttk.Label(tb, text="Width:", font=("TkDefaultFont", 8)).pack(side="left")
        self._rw_var = tk.IntVar(value=row_width)
        sp = ttk.Spinbox(tb, from_=20, to=200, width=5,
                         textvariable=self._rw_var,
                         font=("TkDefaultFont", 8),
                         command=self._on_width_change)
        sp.pack(side="left", padx=2)
        sp.bind("<Return>", lambda e: self._on_width_change())
        ttk.Label(tb, text="bp/row", font=("TkDefaultFont", 8)).pack(side="left")

        # ── legend strip (compact, fixed — does not shrink with window) ─
        leg_f = ttk.Frame(self)
        leg_f.pack(fill="x", pady=(1, 2))
        for ftype, color in UIHelper.COLOR_MAP.items():
            itm = ttk.Frame(leg_f)
            itm.pack(side="left", padx=(0, 4))
            tk.Label(itm, bg=color, width=2, height=1,
                     relief="ridge").pack(side="left")
            ttk.Label(itm, text=ftype,
                      font=("TkDefaultFont", 7)).pack(side="left", padx=(1, 0))

        # ── selection info bar ─────────────────────────────────────────
        info_f = ttk.Frame(self)
        info_f.pack(fill="x")
        self._sel_label = ttk.Label(
            info_f, text="Drag to select  |  Ctrl+C: copy",
            font=("Courier New", 8), foreground="#888888")
        self._sel_label.pack(side="left", padx=4)
        ttk.Button(info_f, text="Copy", width=5,
                   command=self._copy_selection).pack(side="right", padx=2)

        # ── canvas + scrollbars ────────────────────────────────────────
        frm = ttk.Frame(self)
        frm.pack(fill="both", expand=True)
        self.canvas = tk.Canvas(frm, bg="white", highlightthickness=0)
        vsb = ttk.Scrollbar(frm, orient="vertical",   command=self.canvas.yview)
        hsb = ttk.Scrollbar(frm, orient="horizontal", command=self.canvas.xview)
        self.canvas.configure(yscrollcommand=vsb.set, xscrollcommand=hsb.set)
        hsb.pack(side="bottom", fill="x")
        vsb.pack(side="right",  fill="y")
        self.canvas.pack(side="left", fill="both", expand=True)

        # Grab scroll when mouse enters, release when it leaves —
        # this prevents the event from leaking to outer containers.
        self.canvas.bind("<Enter>", self._grab_scroll)
        self.canvas.bind("<Leave>", self._release_scroll)

        self._tip = Tooltip(self.canvas)
        self.canvas.bind("<Motion>",          self._on_motion)
        self.canvas.bind("<Leave>",           lambda e: self._tip.hide())
        self.canvas.bind("<ButtonPress-1>",   self._on_press)
        self.canvas.bind("<B1-Motion>",       self._on_drag)
        self.canvas.bind("<ButtonRelease-1>", self._on_release)
        self.canvas.bind("<Control-c>",       self._copy_selection)
        self.canvas.bind("<Control-C>",       self._copy_selection)

    # ── public API ──────────────────────────────────────────────────────
    def render(self, sequence, features):
        self.sequence = (sequence or "").upper()
        self.features = features or []
        self._draw()

    def show_message(self, msg):
        """Display a plain-text message (for errors / status)."""
        self.canvas.delete("all")
        self._hit_map = []
        self.canvas.create_text(
            10, 10, text=msg, anchor="nw",
            font=("TkDefaultFont", 9), fill="#555555",
            width=600)
        self.canvas.configure(scrollregion=self.canvas.bbox("all") or (0,0,200,50))

    # ── internal ────────────────────────────────────────────────────────
    def _on_width_change(self):
        try:
            v = int(self._rw_var.get())
            self._row_width = max(20, min(200, v))
            
            # Save to settings
            if self.lib:
                self.lib.settings["viewer_row_width"] = self._row_width
                self.lib.save()
                
        except (ValueError, tk.TclError):
            pass
        self._draw()

    def _draw(self):
        c = self.canvas
        c.delete("all")
        self._hit_map = []
        self._row_map = []

        seq = self.sequence
        if not seq:
            return

        rw  = self._row_width
        cw  = self.CW
        ch  = self.CH
        fh  = self.FH
        idx = self.IDX
        total_w = idx + rw * cw + 12

        y = 10
        SEQ_PAD  = 4   # gap between seq row and feature tracks
        LANE_GAP = 3   # gap between feature lanes
        ROW_PAD  = 15  # gap between complete rows

        for row_s in range(0, len(seq), rw):
            row_seq = seq[row_s : row_s + rw]
            row_e   = row_s + len(row_seq)

            # ── 1. Primer tracks (Above) ─────────────────────────────
            row_feats = self._feats_in_range(row_s, row_e)
            primers = [f for f in row_feats if f.get("type") == "Primer"]
            others  = [f for f in row_feats if f.get("type") != "Primer"]

            if primers:
                p_lanes = self._assign_lanes(primers, row_s, row_e)
                for lane in p_lanes:
                    for cs, ce, feat in lane:
                        self._draw_arrow(idx, y, cs, ce, cw, fh, feat)
                    y += fh + LANE_GAP
                y += 2

            # ── 2. Sequence row ──────────────────────────────────────
            y_seq_top = y
            self._row_map.append((row_s, row_e, y_seq_top, y_seq_top + ch))

            # position index
            c.create_text(idx - 6, y + ch // 2,
                          text=str(row_s + 1),
                          anchor="e",
                          font=("Courier New", 8, "italic"),
                          fill="#999999")

            # sequence + background highlights
            for ci, base in enumerate(row_seq):
                abs_pos = row_s + ci
                x = idx + ci * cw
                color = self._color_at(abs_pos)
                if color:
                    c.create_rectangle(x, y, x + cw, y + ch,
                                       fill=color, outline="")
                c.create_text(x + cw // 2, y + ch // 2,
                              text=base,
                              font=("Courier New", 9, "bold") if color else ("Courier New", 9),
                              fill="black")

            # separator line
            c.create_line(idx, y + ch, idx + len(row_seq) * cw, y + ch,
                          fill="#dddddd")

            y += ch + SEQ_PAD

            # ── 3. Other feature tracks (Below) ──────────────────────
            if others:
                o_lanes = self._assign_lanes(others, row_s, row_e)
                for lane in o_lanes:
                    for cs, ce, feat in lane:
                        self._draw_arrow(idx, y, cs, ce, cw, fh, feat)
                    y += fh + LANE_GAP
                y += 2

            y += ROW_PAD

        # update scroll region
        bbox = c.bbox("all")
        if bbox:
            c.configure(scrollregion=(0, 0,
                                      max(bbox[2]+15, total_w),
                                      bbox[3]+15))

        # redraw any active selection on top
        self._draw_selection()

    def _color_at(self, pos):
        # Primers don't highlight the background in SnapGene style (they are on tracks)
        # So we only look for non-primer features for background highlighting.
        best = None
        for f in self.features:
            if f.get("type") == "Primer": continue
            if f.get("start",0) <= pos < f.get("end",0):
                span = f["end"] - f["start"]
                if best is None or span < (best["end"] - best["start"]):
                    best = f
        if best:
            t = best.get("type","Misc")
            return UIHelper.COLOR_MAP.get(t, UIHelper.COLOR_MAP["Misc"])
        return None

    def _feats_in_range(self, rs, re):
        return [f for f in self.features
                if f.get("start",0) < re and f.get("end",0) > rs]

    def _assign_lanes(self, feats, rs, re):
        clipped = []
        for f in feats:
            cs = max(f.get("start",0), rs) - rs
            ce = min(f.get("end",0),  re) - rs
            if ce > cs:
                clipped.append((cs, ce, f))
        clipped.sort(key=lambda x: x[0])
        lanes = []
        for cs, ce, f in clipped:
            placed = False
            for lane in lanes:
                # Add a 1bp buffer between features in the same lane
                if lane[-1][1] + 1 <= cs:
                    lane.append((cs, ce, f))
                    placed = True
                    break
            if not placed:
                lanes.append([(cs, ce, f)])
        return lanes

    def _draw_arrow(self, idx, y, cs, ce, cw, fh, feat):
        x1 = idx + cs * cw
        x2 = idx + ce * cw
        if x2 - x1 < 3: x2 = x1 + 3
        
        ftype  = feat.get("type", "Misc")
        color  = UIHelper.COLOR_MAP.get(ftype, UIHelper.COLOR_MAP["Misc"])
        strand = feat.get("strand", 1)
        label  = feat.get("label", "")
        
        # SnapGene style: Primers have distinctive look
        is_primer = (ftype == "Primer")
        arrow = min(8, max(4, (x2 - x1) // 3))
        my = y + fh // 2

        if strand >= 0:
            pts = [x1, y,  x2-arrow, y,  x2, my,  x2-arrow, y+fh,  x1, y+fh]
        else:
            pts = [x1+arrow, y,  x2, y,  x2, y+fh,  x1+arrow, y+fh,  x1, my]

        outline_color = "#444444" if is_primer else "#777777"
        width = 1.5 if is_primer else 1
        
        self.canvas.create_polygon(pts, fill=color, outline=outline_color, width=width)

        # Label (truncate to fit)
        box_px = x2 - x1 - (arrow if is_primer else 4)
        char_w = self._font_label.measure("A")
        max_chars = max(1, box_px // char_w)
        
        disp = label if len(label) <= max_chars else label[:max_chars-1] + ".."
        self.canvas.create_text((x1+x2)//2, my, text=disp, font=("Courier New", 8), fill="black")

        self._hit_map.append((x1, y, x2, y+fh, feat))

    def _on_motion(self, event):
        if self._dragging:
            return
        cx = self.canvas.canvasx(event.x)
        cy = self.canvas.canvasy(event.y)
        for x1, y1, x2, y2, feat in reversed(self._hit_map):
            if x1 <= cx <= x2 and y1 <= cy <= y2:
                strand = "+" if feat.get("strand", 1) >= 0 else "-"
                tip = "[{}] {} ({}) {}-{}".format(
                    feat.get("type", "Misc"), feat.get("label", ""),
                    strand, feat.get("start", 0)+1, feat.get("end", 0))
                if feat.get("type") == "Primer" and "full_sequence" in feat:
                    tip += f"\nSeq: {feat['full_sequence']}"
                self._tip.show(tip, event.x_root, event.y_root)
                return
        self._tip.hide()

    # ── sequence drag-selection ─────────────────────────────────────────

    def _pos_at(self, cx, cy):
        """Return 0-based sequence position under canvas point (cx, cy), or -1."""
        for row_s, row_e, y_top, y_bot in self._row_map:
            if y_top <= cy <= y_bot:
                col = int((cx - self.IDX) / max(1, self.CW))
                col = max(0, min(col, row_e - row_s - 1))
                return row_s + col
        return -1

    def _draw_selection(self):
        """Draw blue highlight over selected characters (tag='sel')."""
        self.canvas.delete("sel")
        if self._sel_start < 0 or self._sel_end < 0:
            return
        s = min(self._sel_start, self._sel_end)
        e = max(self._sel_start, self._sel_end) + 1
        cw = self.CW
        for row_s, row_e, y_top, y_bot in self._row_map:
            cs = max(s, row_s)
            ce = min(e, row_e)
            if cs >= ce:
                continue
            x1 = self.IDX + (cs - row_s) * cw
            x2 = self.IDX + (ce - row_s) * cw
            self.canvas.create_rectangle(
                x1, y_top, x2, y_bot,
                fill="#3388CC", outline="", stipple="gray25",
                tags="sel")
        self.canvas.tag_raise("sel")
        
        # Update info label with Tm and GC
        sel_seq = self.sequence[s:e]
        bp = len(sel_seq)
        
        tm_str = "N/A"
        gc_str = "N/A"
        if sel_seq and mt and gc_fraction:
            try:
                # Basic Tm calculation using NN model
                tm_val = mt.Tm_NN(sel_seq, Na=50)
                tm_str = f"{tm_val:.1f}°C"
                gc_val = gc_fraction(sel_seq) * 100
                gc_str = f"{gc_val:.1f}%"
            except:
                pass
        
        disp = sel_seq if bp <= 30 else sel_seq[:27] + "..."
        self._sel_label.config(
            text=f"Selection: {s+1}-{e} ({bp} bp) | Tm: {tm_str} | GC: {gc_str} | Seq: {disp}",
            foreground="#1155AA")

    def _on_press(self, event):
        self._tip.hide()
        self._dragging = True
        cx = self.canvas.canvasx(event.x)
        cy = self.canvas.canvasy(event.y)

        # 1. Check if a feature/primer was clicked
        for x1, y1, x2, y2, feat in reversed(self._hit_map):
            if x1 <= cx <= x2 and y1 <= cy <= y2:
                self._sel_start = feat.get("start", 0)
                self._sel_end   = feat.get("end", 1) - 1  # 0-based inclusive
                self._draw_selection()
                self._dragging = False # Clicked feature, don't start dragging
                return

        # 2. Otherwise, start normal drag-selection
        pos = self._pos_at(cx, cy)
        if pos >= 0:
            self._sel_start = pos
            self._sel_end   = pos
            self._draw_selection()
        else:
            self._sel_start = -1
            self._sel_end   = -1
            self.canvas.delete("sel")
            self._sel_label.config(
                text="Drag to select  |  Ctrl+C: copy",
                foreground="#888888")

    def _on_drag(self, event):
        if not self._dragging or self._sel_start < 0:
            return
        cx = self.canvas.canvasx(event.x)
        cy = self.canvas.canvasy(event.y)
        pos = self._pos_at(cx, cy)
        if pos >= 0 and pos != self._sel_end:
            self._sel_end = pos
            self._draw_selection()

    def _on_release(self, event):
        self._dragging = False
        # make canvas focusable so Ctrl+C works
        self.canvas.focus_set()

    def _copy_selection(self, event=None):
        if self._sel_start < 0 or self._sel_end < 0:
            return
        s = min(self._sel_start, self._sel_end)
        e = max(self._sel_start, self._sel_end) + 1
        seq = self.sequence[s:e]
        self.canvas.clipboard_clear()
        self.canvas.clipboard_append(seq)
        self._sel_label.config(
            text="Copied {} bp  (pos {}-{})".format(len(seq), s + 1, e),
            foreground="#117711")

    def _grab_scroll(self, event=None):
        """Mouse entered canvas — take over global scroll so outer frames don't move."""
        self.canvas.bind_all("<MouseWheel>", self._scroll_y)
        self.canvas.bind_all("<Button-4>",   lambda e: self.canvas.yview_scroll(-1, "units"))
        self.canvas.bind_all("<Button-5>",   lambda e: self.canvas.yview_scroll(1,  "units"))

    def _release_scroll(self, event=None):
        """Mouse left canvas — restore default (no bind_all) so outer frames scroll again."""
        self.canvas.unbind_all("<MouseWheel>")
        self.canvas.unbind_all("<Button-4>")
        self.canvas.unbind_all("<Button-5>")

    def _scroll_y(self, event):
        self.canvas.yview_scroll(int(-1 * (event.delta / 120)), "units")


class FindDialog(tk.Toplevel):
    def __init__(self, parent, text_widget, feature_editor=None):
        super().__init__(parent); self.title("Find"); self.geometry("300x80"); self.text_widget = text_widget
        self.feat_ed = feature_editor; self.attributes("-topmost", True)
        ttk.Label(self, text="Pattern:", font=("TkDefaultFont", 8)).pack(pady=2)
        self.entry = ttk.Entry(self, width=30, font=("TkDefaultFont", 8)); self.entry.pack(padx=5); self.entry.focus_set()
        btn_f = ttk.Frame(self); btn_f.pack(pady=5)
        ttk.Button(btn_f, text="Find Next", command=self.find_next).pack(side="left", padx=5)
        self.entry.bind("<Return>", lambda e: self.find_next()); self.last_pos = "1.0"

    def find_next(self):
        p = self.entry.get().strip().upper()
        if not p: return
        self.text_widget.tag_remove("search_hit", "1.0", tk.END)
        self.text_widget.tag_remove("sel", "1.0", tk.END)
        pos = self.text_widget.search(p, self.last_pos, nocase=True, stopindex=tk.END)
        if not pos: pos = self.text_widget.search(p, "1.0", nocase=True, stopindex=tk.END)
        if pos:
            end = f"{pos}+{len(p)}c"
            self.text_widget.tag_add("search_hit", pos, end)
            self.text_widget.tag_add("sel", pos, end)
            self.text_widget.see(pos); self.last_pos = end
            if self.feat_ed: 
                try: 
                    # end 형식을 실제 인덱스(1.15 등)로 정규화한 뒤 추출
                    real_end = self.text_widget.index(end)
                    self.feat_ed.es.delete(0, tk.END); self.feat_ed.es.insert(0, pos.split('.')[1])
                    self.feat_ed.ee.delete(0, tk.END); self.feat_ed.ee.insert(0, real_end.split('.')[1])
                except: pass
