"""Ligation 탭
- Digest fragment / Amplicon / Template 등을 원하는 순서로 배치
- 말단 호환성 자동 체크 (같은 효소 = 호환, Blunt = 호환, 비호환 시 경고)
- 결과 topology를 Linear / Circular 중 선택
- SnapGeneViewer로 결과 미리보기
- 저장 시 category="ligation" 으로 Library에 추가
  → PCR template, HR insert, 추가 Digest source로 활용 가능
"""
import tkinter as tk
from tkinter import ttk, messagebox

from core.models import SequenceItem
from core.restriction import ligate
from gui.gui_components import UIHelper, SearchableCombobox, SnapGeneViewer


class LigationTab:
    def __init__(self, parent, lib, refresh_callback):
        self.parent     = parent
        self.lib        = lib
        self.refresh_cb = refresh_callback
        self.frag_list  = []   # list of dicts: {name, sequence, features, start_enzyme, end_enzyme}
        self.last_result = None
        self.setup_ui()

    # ─────────────────────────────────────────────────────────────
    # UI 구성
    # ─────────────────────────────────────────────────────────────
    def setup_ui(self):
        root = ttk.Frame(self.parent, padding=10)
        root.pack(fill="both", expand=True)

        # ── 1. Fragment 추가 패널 ────────────────────────────────
        add_lf = ttk.LabelFrame(root, text="Fragment 추가", padding=8)
        add_lf.pack(fill="x")

        r0 = ttk.Frame(add_lf); r0.pack(fill="x")
        ttk.Label(r0, text="Source:", width=8).pack(side="left")
        self.src_combo = SearchableCombobox(r0)
        self.src_combo.pack(side="left", fill="x", expand=True, padx=4)
        ttk.Button(r0, text="→ 추가", command=self.add_fragment, width=8).pack(side="left")

        self.src_info = ttk.Label(add_lf, text="(소스 없음)",
                                  font=("TkDefaultFont", 8, "italic"), foreground="gray")
        self.src_info.pack(anchor="w", pady=(2, 0))
        self.src_combo.bind("<<ComboboxSelected>>", lambda e: self._update_src_info())

        # ── 2. Fragment 순서 테이블 ──────────────────────────────
        order_lf = ttk.LabelFrame(root, text="Ligation 순서", padding=4)
        order_lf.pack(fill="x", pady=6)

        tbl_f, self.order_tree = UIHelper.create_scrolled_tree(
            order_lf,
            ("idx", "name", "length", "s_enz", "e_enz"),
            ("#", "Name", "Length (bp)", "5' 말단 효소", "3' 말단 효소"),
        )
        tbl_f.pack(fill="x")
        for col, w in [("idx", 35), ("name", 200), ("length", 90), ("s_enz", 110), ("e_enz", 110)]:
            self.order_tree.column(col, width=w)
        self.order_tree.bind("<<TreeviewSelect>>", lambda e: self._on_frag_select())

        btn_f = ttk.Frame(order_lf); btn_f.pack(fill="x", pady=3)
        ttk.Button(btn_f, text="↑ 위로",   command=self.move_up,      width=8).pack(side="left", padx=2)
        ttk.Button(btn_f, text="↓ 아래로", command=self.move_down,    width=8).pack(side="left", padx=2)
        ttk.Button(btn_f, text="✕ 제거",   command=self.remove_frag,  width=8).pack(side="left", padx=2)
        ttk.Button(btn_f, text="전체 초기화", command=self.clear_all, width=10).pack(side="right", padx=2)

        # ── 3. 호환성 & 옵션 ────────────────────────────────────
        opt_lf = ttk.LabelFrame(root, text="Ligation 옵션 & 호환성", padding=8)
        opt_lf.pack(fill="x", pady=2)

        opt_row = ttk.Frame(opt_lf); opt_row.pack(fill="x")

        ttk.Label(opt_row, text="결과 Topology:", width=16).pack(side="left")
        self.topo_var = tk.StringVar(value="Circular")
        ttk.Radiobutton(opt_row, text="🔵 Circular (환형)",
                        variable=self.topo_var, value="Circular").pack(side="left", padx=6)
        ttk.Radiobutton(opt_row, text="➖ Linear (선형)",
                        variable=self.topo_var, value="Linear").pack(side="left", padx=6)

        ttk.Button(opt_row, text="🔗 Run Ligation",
                   command=self.run_ligation).pack(side="right", padx=4)

        self.compat_label = ttk.Label(opt_lf, text="", font=("TkDefaultFont", 9))
        self.compat_label.pack(anchor="w", pady=(4, 0))

        # ── 4. 결과 뷰어 ────────────────────────────────────────
        res_lf = ttk.LabelFrame(root, text="Ligation 결과 미리보기", padding=4)
        res_lf.pack(fill="both", expand=True, pady=4)

        self.viewer = SnapGeneViewer(res_lf, lib_manager=self.lib)
        self.viewer.pack(fill="both", expand=True)

        # ── 5. 저장 ─────────────────────────────────────────────
        save_f = ttk.Frame(root); save_f.pack(fill="x", pady=4)
        ttk.Label(save_f, text="저장 이름:").pack(side="left")
        self.save_name = ttk.Entry(save_f, width=42)
        self.save_name.pack(side="left", padx=4)
        ttk.Button(save_f, text="💾 Library에 저장",
                   command=self.save_result).pack(side="left", padx=4)

        self.status_lbl = ttk.Label(save_f, text="", font=("TkDefaultFont", 9, "bold"))
        self.status_lbl.pack(side="left", padx=8)

        UIHelper.create_legend(root).pack(side="bottom", fill="x", pady=4)

    # ─────────────────────────────────────────────────────────────
    # 외부 호출: 소스 콤보 갱신
    # ─────────────────────────────────────────────────────────────
    def refresh_sources(self):
        items = []
        for k in sorted(self.lib.digests.keys()):
            v = self.lib.digests[k]
            items.append(f"[Digest] {k}")
        for k in sorted(self.lib.amplicons.keys()):
            items.append(f"[Amplicon] {k}")
        for k in sorted(self.lib.ligations.keys()):
            items.append(f"[Ligation] {k}")
        for k in sorted(self.lib.templates.keys()):
            items.append(f"[Template] {k}")
        self.src_combo['values'] = items

    # ─────────────────────────────────────────────────────────────
    # 내부 헬퍼
    # ─────────────────────────────────────────────────────────────
    def _get_src_obj(self):
        """콤보박스에서 선택된 항목의 SequenceItem 반환."""
        s = self.src_combo.get()
        if not s:
            return None
        tag, name = (s.split("] ", 1) + [""])[:2] if "] " in s else ("", s)
        tag = tag.lstrip("[")
        store = {
            "Digest":   self.lib.digests,
            "Amplicon": self.lib.amplicons,
            "Ligation": self.lib.ligations,
            "Template": self.lib.templates,
        }.get(tag, {})
        return store.get(name)

    def _update_src_info(self):
        obj = self._get_src_obj()
        if obj:
            s_enz = getattr(obj, "enzymes", "") or "?"
            topo  = getattr(obj, "topology", "")
            self.src_info.config(
                text=f"길이: {len(obj.sequence)} bp  |  Topology: {topo}"
                     f"  |  효소: {s_enz if s_enz != '?' else '정보 없음'}",
                foreground="black"
            )
        else:
            self.src_info.config(text="(소스 없음)", foreground="gray")

    def _frag_from_obj(self, obj, display_name):
        """SequenceItem → ligation용 dict로 변환."""
        # enzymes 필드: "EcoRI,SalI" 형태 → 첫 효소가 5' 말단, 마지막이 3' 말단
        enzymes = getattr(obj, "enzymes", "") or ""
        enz_list = [e.strip() for e in enzymes.split(",") if e.strip()]
        s_enz = enz_list[0]  if enz_list else "Blunt"
        e_enz = enz_list[-1] if enz_list else "Blunt"
        return {
            "name":         display_name,
            "sequence":     obj.sequence,
            "features":     list(obj.features) if obj.features else [],
            "start_enzyme": s_enz,
            "end_enzyme":   e_enz,
            "topology":     getattr(obj, "topology", "Linear"),
        }

    def _rebuild_table(self):
        self.order_tree.delete(*self.order_tree.get_children())
        for i, f in enumerate(self.frag_list):
            self.order_tree.insert("", "end", iid=str(i), values=(
                i + 1,
                f["name"],
                len(f["sequence"]),
                f["start_enzyme"],
                f["end_enzyme"],
            ))
        self._update_compat_display()

    def _update_compat_display(self):
        """Fragment 목록이 바뀔 때마다 말단 호환성 요약을 갱신."""
        n = len(self.frag_list)
        if n < 2:
            self.compat_label.config(text="Fragment를 2개 이상 추가하세요.", foreground="gray")
            return

        issues = []
        check_pairs = list(range(n - 1))
        if self.topo_var.get() == "Circular":
            check_pairs.append(-1)   # 마지막 ↔ 첫 번째

        from core.restriction import _ends_compatible
        for i in check_pairs:
            left  = self.frag_list[i]
            right = self.frag_list[(i + 1) % n]
            le, re = left["end_enzyme"], right["start_enzyme"]
            if not _ends_compatible(le, re):
                issues.append(f"{left['name']}(3':{le}) ↔ {right['name']}(5':{re})")

        if issues:
            self.compat_label.config(
                text="⚠️ 비호환 접합부: " + " | ".join(issues),
                foreground="#CC4400"
            )
        else:
            self.compat_label.config(
                text=f"✅ 모든 접합부 호환 ({n - 1 + (1 if self.topo_var.get() == 'Circular' else 0)}개 junction)",
                foreground="#006600"
            )

    def _on_frag_select(self):
        pass  # 향후 확장용 (선택된 fragment preview 등)

    def _suggest_name(self):
        if not self.frag_list:
            return ""
        names = [f["name"].split("]")[-1].strip() for f in self.frag_list]
        topo  = self.topo_var.get()
        return f"Ligation ({' + '.join(names)}) [{topo}]"

    # ─────────────────────────────────────────────────────────────
    # 버튼 핸들러
    # ─────────────────────────────────────────────────────────────
    def add_fragment(self):
        s = self.src_combo.get()
        if not s:
            messagebox.showwarning("경고", "소스를 선택하세요.")
            return
        obj = self._get_src_obj()
        if obj is None:
            messagebox.showwarning("경고", f"'{s}' 를 찾을 수 없습니다.")
            return
        frag = self._frag_from_obj(obj, s)
        self.frag_list.append(frag)
        self._rebuild_table()
        # 자동 이름 제안
        self.save_name.delete(0, tk.END)
        self.save_name.insert(0, self._suggest_name())

    def move_up(self):
        sel = self.order_tree.selection()
        if not sel:
            return
        i = int(sel[0])
        if i == 0:
            return
        self.frag_list[i], self.frag_list[i - 1] = self.frag_list[i - 1], self.frag_list[i]
        self._rebuild_table()
        self.order_tree.selection_set(str(i - 1))

    def move_down(self):
        sel = self.order_tree.selection()
        if not sel:
            return
        i = int(sel[0])
        if i >= len(self.frag_list) - 1:
            return
        self.frag_list[i], self.frag_list[i + 1] = self.frag_list[i + 1], self.frag_list[i]
        self._rebuild_table()
        self.order_tree.selection_set(str(i + 1))

    def remove_frag(self):
        sel = self.order_tree.selection()
        if not sel:
            return
        i = int(sel[0])
        del self.frag_list[i]
        self._rebuild_table()

    def clear_all(self):
        self.frag_list.clear()
        self._rebuild_table()
        self.viewer.show_message("Fragment를 추가하고 Run Ligation을 실행하세요.")
        self.last_result = None
        self.status_lbl.config(text="")

    def run_ligation(self):
        if len(self.frag_list) < 2:
            messagebox.showwarning("경고", "Fragment를 2개 이상 추가하세요.")
            return

        topo = self.topo_var.get()

        # 비호환 junction 경고 (진행은 허용)
        from core.restriction import _ends_compatible
        n = len(self.frag_list)
        bad = []
        pairs = list(range(n - 1))
        if topo == "Circular":
            pairs.append(-1)
        for i in pairs:
            left  = self.frag_list[i]
            right = self.frag_list[(i + 1) % n]
            if not _ends_compatible(left["end_enzyme"], right["start_enzyme"]):
                bad.append(f"{left['name']} ↔ {right['name']}")
        if bad:
            proceed = messagebox.askyesno(
                "호환성 경고",
                "비호환 말단이 있습니다:\n" + "\n".join(bad)
                + "\n\n그래도 Ligation을 진행하시겠습니까?\n"
                  "(강제 결합 — 실제 실험에서는 효율이 낮을 수 있습니다)"
            )
            if not proceed:
                return

        result = ligate(self.frag_list, topology=topo)
        if result is None:
            messagebox.showerror("오류", "Ligation에 실패했습니다.")
            return

        self.last_result = result

        # Junction을 feature로 추가하여 시각화
        viz_feats = list(result["features"])
        for j in result["junctions"]:
            pos = j["pos"]
            if pos == 0 and topo == "Circular":
                continue   # 원형의 position-0 junction은 표시 생략
            viz_feats.append({
                "label":  f"Junction ({j['left']}|{j['right']})",
                "start":  max(0, pos - 1),
                "end":    min(len(result["sequence"]), pos + 1),
                "type":   "Misc",
                "strand": 1,
            })

        UIHelper.render_gene_viewer(self.viewer, result["sequence"], viz_feats,
                                    lib_manager=self.lib)

        self.status_lbl.config(
            text=f"결과: {len(result['sequence'])} bp  [{topo}]",
            foreground="darkblue"
        )
        # 이름 제안 업데이트
        self.save_name.delete(0, tk.END)
        self.save_name.insert(0, self._suggest_name())

    def save_result(self):
        if self.last_result is None:
            messagebox.showwarning("경고", "먼저 Run Ligation을 실행하세요.")
            return
        name = self.save_name.get().strip()
        if not name:
            messagebox.showwarning("경고", "저장 이름을 입력하세요.")
            return

        # 이름 중복 회피
        original = name
        ctr = 2
        while name in self.lib.ligations:
            name = f"{original} ({ctr})"
            ctr += 1

        src_names = ", ".join(f["name"] for f in self.frag_list)
        item = SequenceItem(
            name=name,
            sequence=self.last_result["sequence"],
            category="ligation",
            topology=self.last_result["topology"],
            features=self.last_result["features"],
            template_name=src_names,
            enzymes=",".join(
                sorted(set(
                    j["left"] for j in self.last_result["junctions"]
                ) | set(
                    j["right"] for j in self.last_result["junctions"]
                ))
            ),
        )
        self.lib.ligations[name] = item
        self.lib.save()
        self.refresh_cb()

        self.status_lbl.config(
            text=f"✅ '{name}' 저장 완료!",
            foreground="#006600"
        )
        messagebox.showinfo(
            "저장 완료",
            f"'{name}'\n길이: {len(item.sequence)} bp  [{item.topology}]\n\n"
            "Library > Ligation Product 탭에서 확인하세요."
        )
