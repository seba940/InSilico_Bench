import tkinter as tk
from tkinter import ttk
from models import GlobalFeature, SequenceItem, Marker, EpitopeTag
from manager import LibraryManager
from gui_components import UIHelper

def test_cassette_viz():
    root = tk.Tk()
    root.title("Cassette Visualization Test")
    root.geometry("1000x600")
    
    lib = LibraryManager("test_data.json")
    
    # 1. 테스트용 전역 피처 등록
    lib.global_features["pTEF1"] = GlobalFeature("pTEF1", "GATCC" * 10, "Promoter")
    lib.global_features["URA3"] = GlobalFeature("URA3", "ATGCG" * 20, "Marker")
    lib.global_features["tCYC1"] = GlobalFeature("tCYC1", "TTAAA" * 10, "Terminator")
    lib.global_features["3xHA"] = GlobalFeature("3xHA", "ATGCGATGCG", "Tag") # 중첩 테스트용
    
    # 2. 카세트 서열 생성
    cassette_seq = ("GATCC" * 10) + ("ATGCG" * 20) + ("TTAAA" * 10)
    
    # 3. UI 구성
    main_f = ttk.Frame(root, padding=20)
    main_f.pack(expand=True, fill="both")
    
    ttk.Label(main_f, text="Cassette: pTEF1 - URA3 - tCYC1", font=("TkDefaultFont", 12, "bold")).pack(pady=10)
    
    viewer = tk.Text(main_f, height=10, wrap="char", font=("Courier", 12))
    viewer.pack(expand=True, fill="both")
    
    # 범례 추가
    UIHelper.create_legend(main_f).pack(fill="x", pady=10)
    
    # 4. 렌더링 실행
    UIHelper.render_annotations(viewer, cassette_seq, [], lib_manager=lib)
    
    root.mainloop()

if __name__ == "__main__":
    test_cassette_viz()
