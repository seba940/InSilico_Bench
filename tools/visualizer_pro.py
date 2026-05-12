import streamlit as st
from Bio.SeqFeature import SeqFeature, FeatureLocation
from dataclasses import dataclass
from typing import List, Dict

# 1. 종류별 색상 맵핑 (Color Mapping Strategy)
COLOR_MAP = {
    "Promoter": "#B2FFB2",      # 연초록
    "CDS": "#B2D8FF",           # 연파랑
    "Marker": "#FFB2B2",        # 연빨강 (Selection Marker)
    "Tag": "#E6B2FF",           # 연보라 (Epitope Tag)
    "Terminator": "#FFFFB2",    # 연노랑
    "Homology Arm": "#E0E0E0",  # 회색
    "Misc": "#F5DEB3"           # 연갈색
}

@dataclass
class GlobalFeature:
    label: str
    start: int
    end: int
    type: str
    strand: int = 1

def generate_html_visualization(sequence: str, features: List[GlobalFeature]):
    """
    서열 바탕색 렌더링 (Sequence Background Highlighting)
    지능형 중첩 충돌 방지 로직 포함 (좁은 범위 우선)
    """
    seq_len = len(sequence)
    bp_colors = [None] * seq_len
    bp_labels = [None] * seq_len
    
    # 2. 중첩 처리 로직: 긴 Feature부터 먼저 칠하고, 짧은 Feature를 나중에 덮어씌움
    sorted_feats = sorted(features, key=lambda x: (x.end - x.start), reverse=True)
    
    for f in sorted_feats:
        color = COLOR_MAP.get(f.type, COLOR_MAP["Misc"])
        for i in range(max(0, f.start), min(seq_len, f.end)):
            bp_colors[i] = color
            bp_labels[i] = f.label

    # HTML 생성
    html_res = '<div style="font-family: monospace; white-space: pre-wrap; line-height: 1.8; font-size: 14px; background-color: #f8f9fa; padding: 20px; border-radius: 10px; border: 1px solid #ddd; word-break: break-all;">'
    
    current_color = None
    for i, base in enumerate(sequence):
        color = bp_colors[i]
        label = bp_labels[i]
        
        if color != current_color:
            if current_color is not None:
                html_res += '</span>'
            if color is not None:
                html_res += f'<span title="{label}" style="background-color: {color}; border-radius: 2px; padding: 2px 0; border: 0.5px solid rgba(0,0,0,0.1);">'
            current_color = color
        
        html_res += base
        
    if current_color is not None:
        html_res += '</span>'
    
    html_res += '</div>'
    return html_res

def render_legend():
    """3. 주석 색상 범례 (Annotation Legend) 패널"""
    st.write("### 🎨 Annotation Legend")
    # 범례를 한 줄에 표시
    legend_html = '<div style="display: flex; flex-wrap: wrap; gap: 15px; background: white; padding: 10px; border-radius: 5px; border: 1px solid #eee;">'
    for f_type, color in COLOR_MAP.items():
        legend_html += f"""
            <div style="display: flex; align-items: center; gap: 6px;">
                <div style="width: 18px; height: 18px; background-color: {color}; border-radius: 3px; border: 1px solid #999;"></div>
                <span style="font-size: 13px; font-weight: 500;">{f_type}</span>
            </div>
        """
    legend_html += '</div>'
    st.markdown(legend_html, unsafe_allow_html=True)

# 4. 통합 테스트 (Integrated Test)
if __name__ == "__main__":
    st.set_page_config(layout="wide", page_title="DNA Pro Visualizer")
    st.title("🧬 DNA Pro Visualizer")
    
    # 샘플 서열 및 데이터 생성
    ptef1_seq = "GATCC" * 10 # 50bp
    ura3_seq = "ATGCG" * 20 # 100bp
    tcyc1_seq = "TTAAA" * 10 # 50bp
    full_sequence = ptef1_seq + ura3_seq + tcyc1_seq
    
    test_features = [
        GlobalFeature("pTEF1", 0, 50, "Promoter"),
        GlobalFeature("URA3", 50, 150, "Marker"),
        GlobalFeature("tCYC1", 150, 200, "Terminator"),
        GlobalFeature("3xHA Tag", 80, 110, "Tag") # 중첩 테스트
    ]
    
    st.subheader("Sequence Annotation Map")
    html_map = generate_html_visualization(full_sequence, test_features)
    st.markdown(html_map, unsafe_allow_html=True)
    
    st.markdown("<br>", unsafe_allow_html=True)
    render_legend()
    
    st.success("✅ SnapGene-style visualization is active!")
