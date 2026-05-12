# Yeast Lab Simulator (Pro v14)

Yeast Lab Simulator는 효모(Saccharomyces cerevisiae) 실험실에서 수행하는 주요 분자생물학 실험을 디지털로 시뮬레이션하고 관리할 수 있는 전문 도구입니다. 복잡한 PCR 설계, 상동 재조합(HR), 제한효소 처리 과정을 직관적인 UI와 SnapGene 스타일의 서열 뷰어를 통해 미리 확인하고 최적화할 수 있습니다.

## 주요 기능

### 1. 🧬 PCR Execution (PCR 시뮬레이션)
*   **주형 및 프라이머 결합 분석**: 선택한 주형(Template)과 프라이머 쌍의 결합 위치 및 방향성을 자동 탐색합니다.
*   **SnapGene 스타일 시각화**: 주형 서열 위에 프라이머가 결합하는 위치와 방향을 화살표 트랙으로 표시합니다.
*   **결과 예측**: PCR 산물(Amplicon)의 예상 길이, 어닐링 온도(Ta), 신장 시간(Extension Time)을 자동으로 계산합니다.
*   **라이브러리 저장**: 생성된 Amplicon을 라이브러리에 저장하여 HR 시뮬레이션의 Insert로 즉시 활용 가능합니다.

### 2. 🔄 HR Simulation (상동 재조합 시뮬레이션)
*   **유전자 교체/삽입**: Target Genome과 Insert(Amplicon/Digest fragment) 간의 상동 재조합 결과를 예측합니다.
*   **자동 말단 탐색**: 양 끝단의 상동 영역(Homology Arm)을 자동으로 찾아 재조합된 최종 서열을 생성합니다.
*   **Feature 상속**: 재조합 과정에서 주석(Feature), 마커, 태그 정보가 유실되지 않고 최종 서열에 유지됩니다.

### 3. ✂️ Restriction Digest (제한효소 처리)
*   **REBASE 데이터베이스**: 수백 종의 제한효소 정보를 바탕으로 서열 내 Cut Site를 검색합니다.
*   **다중 효소 지원**: 여러 효소를 동시에 처리했을 때 생성되는 모든 단편(Fragments)을 목록화합니다.
*   **단편 미리보기**: 각 단편의 서열과 포함된 Feature를 SnapGene 뷰어로 즉시 확인할 수 있습니다.

### 4. 🧬 Gene Lookup (유전자 검색 및 추출)
*   **SGD 연동**: 유전자 이름이나 별칭으로 효모 유전체 데이터베이스를 검색합니다.
*   **상/하류(Flanking) 추출**: 유전자 본체뿐만 아니라 원하는 길이의 Upstream/Downstream 영역을 포함하여 서열을 추출합니다.
*   **SGD 페이지 바로가기**: 선택한 유전자의 상세 정보를 웹 브라우저에서 즉시 확인할 수 있습니다.

### 5. 📁 Library Management (체계적인 라이브러리 관리)
*   **통합 관리**: 프라이머, 주형, 마커, 태그, 종(Species), 주석 참조 파일을 한곳에서 관리합니다.
*   **프라이머 검증**: 단일 프라이머의 이차 구조(Hairpin, Dimer) 분석 및 Tm 값을 계산합니다.
*   **가변 레이아웃**: 프라이머 목록 등 주요 테이블의 너비를 자유롭게 조절할 수 있습니다.

---

## 설치 및 실행 방법 (개발자용)

이 프로그램은 Python 3.10 이상 환경에서 구동됩니다.

### 1. 필수 라이브러리 설치
터미널 또는 명령 프롬프트에서 다음 명령어를 실행하여 필요한 라이브러리를 설치합니다.
```bash
pip install biopython pandas primer3-py openpyxl
```

### 2. 프로그램 실행
프로젝트 루트 폴더에서 다음 파일을 실행합니다.
```bash
python main.py
```

---

## 기술 스택
*   **Language**: Python 3.10+
*   **GUI Framework**: Tkinter (ttk)
*   **Bioinformatics**: Biopython, Primer3-py
*   **Data Management**: Pandas, JSON

## 라이선스
본 소프트웨어는 연구 및 교육 목적으로 제작되었습니다.
