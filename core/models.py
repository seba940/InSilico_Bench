from dataclasses import dataclass, field
from typing import List, Dict

@dataclass
class Primer:
    name: str = ""
    full_sequence: str = ""
    binding_len: int = 20
    p_type: str = "Fwd"
    usage: str = "PCR"
    tm: float = 0.0        # Binding site Tm
    tm_total: float = 0.0  # Full sequence Tm

    def __post_init__(self):
        try:
            try:
                from Bio.SeqUtils.MeltingTemp import Tm_NN
            except ImportError:
                from Bio.SeqUtils import Tm_NN
            
            # Binding site Tm
            if not self.tm:
                bind_seq = self.full_sequence[-self.binding_len:]
                self.tm = round(Tm_NN(bind_seq, Na=50, dnac1=250, dnac2=250), 1)
            
            # Total sequence Tm
            if not self.tm_total:
                self.tm_total = round(Tm_NN(self.full_sequence, Na=50, dnac1=250, dnac2=250), 1)
        except Exception:
            self.tm = self.tm or 0.0
            self.tm_total = self.tm_total or 0.0

    @property
    def binding(self):
        return self.full_sequence[-self.binding_len:]

    @property
    def overhang(self):
        return self.full_sequence[:-self.binding_len]


@dataclass
class SequenceItem:
    name: str = ""
    sequence: str = ""
    category: str = "template"
    topology: str = "Linear"
    species: str = ""
    features: List[Dict] = field(default_factory=list)
    template_name: str = ""
    marker: str = ""
    tag: str = ""
    enzymes: str = ""  # digest 카테고리에서 사용된 효소 (콤마 구분)
    kind: str = ""     # template 한정: "Genome" 또는 "Plasmid"


@dataclass
class Marker:
    name: str = ""
    sequence: str = ""
    features: List[Dict] = field(default_factory=list)


@dataclass
class EpitopeTag:
    name: str = ""
    sequence: str = ""
    features: List[Dict] = field(default_factory=list)


@dataclass
class Species:
    name: str = ""
    description: str = ""


@dataclass
class AnnotationRef:
    name: str = ""
    features: List[Dict] = field(default_factory=list)
    source_file: str = ""
    species: str = "Generic"
    provider: str = "NCBI"


@dataclass
class GlobalFeature:
    label: str = ""
    sequence: str = ""
    type: str = "Misc"
    strand: int = 1
