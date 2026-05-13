from Bio.Seq import Seq
from Bio.SeqUtils import MeltingTemp as mt

class PrimerManager:
    def __init__(self):
        # 프라이머 데이터를 메모리에 저장할 딕셔너리
        self.primers = {}

    def calculate_tm(self, sequence):
        """
        주어진 서열의 Tm 값을 Nearest-Neighbor 방식으로 계산합니다.
        """
        my_seq = Seq(sequence.upper())
        try:
            # SnapGene과 유사한 결과를 얻기 위해 Tm_NN 사용
            # 필요에 따라 Na+(염) 농도, DNA 농도 파라미터를 세밀하게 조정할 수 있습니다.
            tm_value = mt.Tm_NN(my_seq, Na=50) # 일반적인 50mM Na+ 조건 가정
            return round(tm_value, 2)
        except Exception as e:
            print(f"Tm 계산 중 오류 발생: {e}")
            return None

    def add_primer(self, name, sequence, description=""):
        """
        프라이머 이름, 서열, 설명(선택)을 입력받아 저장합니다.
        """
        sequence = sequence.upper()
        tm = self.calculate_tm(sequence)
        
        self.primers[name] = {
            "sequence": sequence,
            "tm": tm,
            "length": len(sequence),
            "description": description
        }
        print(f"[등록 완료] {name} | 길이: {len(sequence)}bp | Tm: {tm}°C")

    def get_primer(self, name):
        """
        저장된 프라이머 정보를 불러옵니다.
        """
        if name in self.primers:
            return self.primers[name]
        else:
            print(f"[오류] '{name}' 프라이머를 찾을 수 없습니다.")
            return None

    def show_all(self):
        """
        현재 등록된 모든 프라이머 목록을 출력합니다.
        """
        print(f"\n=== 저장된 프라이머 목록 (총 {len(self.primers)}개) ===")
        for name, info in self.primers.items():
            print(f" - {name}: {info['sequence']} (Tm: {info['tm']}°C)")
        print("=========================================\n")

# --- 테스트 실행 코드 ---
if __name__ == "__main__":
    pm = PrimerManager()
    
    # 예시 프라이머 등록
    pm.add_primer("Fwd_KO_Primer", "CGTACGTAGCTAGCTAGCTAGCTAGCATCGATCGATCGA", "Target Gene KO Forward")
    pm.add_primer("Rev_Check", "GATCGATCGATCGATCGATC", "Genotyping Reverse")
    
    # 목록 확인
    pm.show_all()

    # 창이 바로 닫히는 것을 방지하는 코드 (이 줄을 추가하세요!)
    input("결과를 확인하셨으면 엔터 키를 눌러 종료하세요...")