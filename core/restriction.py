"""제한효소 처리 헬퍼 모듈
- Biopython의 Bio.Restriction (REBASE) 전체 enzyme 목록 사용
- cut site 검색, fragment 생성, feature 상속 로직 제공
"""
from Bio.Seq import Seq
from Bio.Restriction import AllEnzymes, RestrictionBatch, Restriction


def get_all_enzyme_names():
    """REBASE 전체 enzyme 이름을 알파벳 순으로 반환."""
    try:
        return sorted([str(e) for e in AllEnzymes])
    except Exception:
        return []


def get_enzyme(name):
    """이름으로 enzyme 객체 반환 (없으면 None)."""
    try:
        return getattr(Restriction, name, None)
    except Exception:
        return None


def enzyme_info(name):
    """enzyme의 recognition site, 절단부 정보 등 텍스트 요약 반환."""
    e = get_enzyme(name)
    if e is None: return f"{name}: (not found)"
    try:
        return f"{name}: site={e.site}, length={e.size}, ovhg={e.ovhg}"
    except Exception:
        return name


def find_cut_sites(sequence, topology, enzyme_names):
    """선택된 enzyme들로 sequence에서 cut site 검색.
    반환: [{"pos": int(0-based cut position on top strand),
           "enzyme": str,
           "site": str(recognition site)}, ...]  pos 오름차순 정렬
    """
    if not sequence or not enzyme_names: return []
    seq = Seq(sequence)
    is_linear = (topology.lower() != "circular")

    # 유효한 enzyme만 모음
    valid_names = [n for n in enzyme_names if get_enzyme(n) is not None]
    if not valid_names: return []

    sites = []
    try:
        rb = RestrictionBatch(valid_names)
        result = rb.search(seq, linear=is_linear)
        # result: {EnzymeClass: [1-based positions on top strand], ...}
        for enz, positions in result.items():
            ename = str(enz)
            esite = getattr(enz, "site", "")
            for p in positions:
                # Biopython은 cut 직후 base의 1-based 위치를 반환
                # → top strand에서의 0-based cut position = p - 1
                cut_pos = (p - 1) % len(sequence)
                sites.append({"pos": cut_pos, "enzyme": ename, "site": esite})
    except Exception as ex:
        # 일부 enzyme이 잘못됐을 수 있으니 개별로 폴백
        for n in valid_names:
            e = get_enzyme(n)
            if not e: continue
            try:
                positions = e.search(seq, linear=is_linear)
                esite = getattr(e, "site", "")
                for p in positions:
                    sites.append({"pos": (p - 1) % len(sequence),
                                  "enzyme": n, "site": esite})
            except Exception:
                continue

    sites.sort(key=lambda x: x["pos"])
    return sites


def digest(sequence, topology, sites):
    """cut site 목록으로 sequence를 잘라 fragment list 반환.
    각 fragment dict: {
      "sequence": str,
      "cut_start": int,   # 원본 sequence 상의 시작 cut 위치 (0-based)
      "cut_end": int,     # 원본 sequence 상의 끝 cut 위치 (0-based)
      "wrapped": bool,    # 원본을 wrap-around 했는지 여부
      "start_enzyme": str,
      "end_enzyme": str,
    }
    """
    if not sites: return []
    seq_len = len(sequence)
    is_circular = (topology.lower() == "circular")

    # 위치별 enzyme 이름 매핑 (같은 위치에 여러 enzyme 가능)
    enz_at = {}
    for s in sites:
        enz_at.setdefault(s["pos"], []).append(s["enzyme"])
    cuts = sorted(enz_at.keys())

    fragments = []
    if is_circular:
        # circular: cuts[i] ~ cuts[i+1]까지의 fragment (마지막은 wrap)
        n = len(cuts)
        for i in range(n):
            cs = cuts[i]
            ce = cuts[(i + 1) % n]
            if ce > cs:
                seq_frag = sequence[cs:ce]
                wrapped = False
            else:
                seq_frag = sequence[cs:] + sequence[:ce]
                wrapped = True
            fragments.append({
                "sequence": seq_frag,
                "cut_start": cs, "cut_end": ce, "wrapped": wrapped,
                "start_enzyme": ",".join(sorted(set(enz_at[cs]))),
                "end_enzyme": ",".join(sorted(set(enz_at[ce]))),
            })
        # 단일 cut의 경우 위 루프는 wrap-around 한 fragment 1개를 반환
    else:
        # linear: 0~cuts[0], cuts[i]~cuts[i+1], cuts[-1]~end
        prev = 0
        prev_enz = "5'-end"
        for p in cuts:
            fragments.append({
                "sequence": sequence[prev:p],
                "cut_start": prev, "cut_end": p, "wrapped": False,
                "start_enzyme": prev_enz,
                "end_enzyme": ",".join(sorted(set(enz_at[p]))),
            })
            prev = p
            prev_enz = ",".join(sorted(set(enz_at[p])))
        fragments.append({
            "sequence": sequence[prev:seq_len],
            "cut_start": prev, "cut_end": seq_len, "wrapped": False,
            "start_enzyme": prev_enz, "end_enzyme": "3'-end",
        })

    return fragments


def inherit_features(features, fragment, source_topology, source_len):
    """원본 feature를 fragment 좌표계로 변환. cut을 가로지르는 feature는 분할."""
    if not features: return []
    cs, ce = fragment["cut_start"], fragment["cut_end"]
    new_feats = []
    is_wrapped = fragment.get("wrapped", False)

    if is_wrapped:
        # fragment는 [cs, source_len) + [0, ce) 의 결합
        # 좌표 변환: orig_idx in [cs, source_len) → frag_idx = orig_idx - cs
        #           orig_idx in [0, ce) → frag_idx = source_len - cs + orig_idx
        offset_part2 = source_len - cs
        for f in features:
            fs, fe = f["start"], f["end"]
            in_part1 = (fs >= cs and fe <= source_len)
            in_part2 = (fs >= 0 and fe <= ce)
            if in_part1:
                nf = f.copy(); nf["start"] = fs - cs; nf["end"] = fe - cs
                new_feats.append(nf)
            elif in_part2:
                nf = f.copy(); nf["start"] = fs + offset_part2; nf["end"] = fe + offset_part2
                new_feats.append(nf)
            # cut을 가로지르는 feature는 양쪽으로 분할 시도
            elif fs < cs < fe and fe <= source_len:
                # part1 쪽 일부만 살림
                nf = f.copy(); nf["label"] = f.get("label", "") + "*"
                nf["start"] = 0; nf["end"] = fe - cs
                if nf["end"] - nf["start"] > 0: new_feats.append(nf)
            elif fs < ce and fe > ce and fs >= 0:
                # part2 쪽 일부만 살림
                nf = f.copy(); nf["label"] = f.get("label", "") + "*"
                nf["start"] = fs + offset_part2
                nf["end"] = ce + offset_part2
                if nf["end"] - nf["start"] > 0: new_feats.append(nf)
    else:
        # 단순 [cs, ce) 구간
        for f in features:
            fs, fe = f["start"], f["end"]
            if fs >= cs and fe <= ce:
                nf = f.copy(); nf["start"] = fs - cs; nf["end"] = fe - cs
                new_feats.append(nf)
            elif fs < cs < fe <= ce:
                # cs를 가로지르는 feature → 우측 일부만
                nf = f.copy(); nf["label"] = f.get("label", "") + "*"
                nf["start"] = 0; nf["end"] = fe - cs
                if nf["end"] - nf["start"] > 0: new_feats.append(nf)
            elif cs <= fs < ce < fe:
                # ce를 가로지르는 feature → 좌측 일부만
                nf = f.copy(); nf["label"] = f.get("label", "") + "*"
                nf["start"] = fs - cs; nf["end"] = ce - cs
                if nf["end"] - nf["start"] > 0: new_feats.append(nf)

    return new_feats


def annotate_cut_sites_on_source(sequence, sites):
    """원본 sequence에 cut site들을 feature로 표시 (시각화용).
    각 cut을 길이 1짜리 'Misc' 타입 feature로 등록.
    """
    feats = []
    for s in sites:
        p = s["pos"]
        feats.append({
            "label": f"{s['enzyme']} cut",
            "start": p,
            "end": min(p + 1, len(sequence)),
            "type": "Misc",
            "strand": 1,
        })
    return feats
