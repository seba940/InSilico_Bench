from Bio.Seq import Seq
from Bio import SeqIO
import os
import re
import urllib.parse

# IUPAC degenerate base regex mapping
_IUPAC_REGEX = {
    'R': '[AG]', 'Y': '[CT]', 'S': '[GC]', 'W': '[AT]',
    'K': '[GT]', 'M': '[AC]', 'B': '[CGT]', 'D': '[AGT]',
    'H': '[ACT]', 'V': '[ACG]', 'N': '[ACGT]',
}
_VALID_DNA = set('ACGTNRYSWKMBDHV')

# Extension time in seconds per 1000 bp
POLYMERASE_RATES = {
    "Taq":     60,
    "Q5":      30,
    "Phusion": 30,
    "KOD":     20,
    "Pfu":     120,
}


def validate_dna_sequence(seq):
    """Return list of invalid characters. Empty list means valid."""
    return [c for c in seq.upper() if c not in _VALID_DNA]


def get_rc(sequence):
    return str(Seq(sequence).reverse_complement())


def _iupac_to_regex(seq):
    pattern = ''
    for c in seq.upper():
        if c in _IUPAC_REGEX:
            pattern += _IUPAC_REGEX[c]
        else:
            pattern += re.escape(c)
    return pattern


def _has_degen(seq):
    return any(c in 'RYSWKMBDHVN' for c in seq.upper())


def smart_find_binding(template, primer_binding, p_type="Fwd", topology="Linear"):
    temp_fwd = template.upper()
    pb = primer_binding.upper()
    degen = _has_degen(pb)

    if p_type == "Fwd":
        search_temp = temp_fwd + temp_fwd[:len(pb)] if topology == "Circular" else temp_fwd
        for length in range(len(pb), 9, -1):
            sub = pb[-length:]
            if degen:
                m = re.search(_iupac_to_regex(sub), search_temp)
                if m:
                    return (m.start() % len(temp_fwd)), length
            else:
                idx = search_temp.find(sub)
                if idx != -1:
                    return (idx % len(temp_fwd)), length
    else:
        temp_rc = get_rc(temp_fwd)
        orig_len = len(temp_fwd)
        search_temp = temp_rc + temp_rc[:len(pb)] if topology == "Circular" else temp_rc
        for length in range(len(pb), 9, -1):
            sub = pb[-length:]
            if degen:
                m = re.search(_iupac_to_regex(sub), search_temp)
                if m:
                    actual_idx_rc = m.start() % orig_len
                    idx_fwd = (orig_len - actual_idx_rc - length) % orig_len
                    return idx_fwd, length
            else:
                idx_rc = search_temp.find(sub)
                if idx_rc != -1:
                    actual_idx_rc = idx_rc % orig_len
                    idx_fwd = (orig_len - actual_idx_rc - length) % orig_len
                    return idx_fwd, length
    return -1, 0


def find_primer_bindings(sequence, fwd_primers, rev_primers):
    """Search every primer's binding region in *sequence* (both strands).

    Returns a list of feature dicts:
      {"label": name, "start": int, "end": int,
       "type": "Primer", "strand": +1 or -1}

    + strand match  → strand = +1  (right-pointing arrow in viewer)
    - strand match  → strand = -1  (left-pointing arrow in viewer)

    Only the binding portion (3' end, ``binding_len`` bp) is searched,
    not the full primer including overhang.
    """
    results = []
    seq_up = sequence.upper()
    seen = set()  # (start, end, strand) dedup

    all_primers = []
    if isinstance(fwd_primers, dict):
        all_primers.extend(list(fwd_primers.values()))
    elif isinstance(fwd_primers, (list, tuple)):
        all_primers.extend(fwd_primers)
        
    if isinstance(rev_primers, dict):
        all_primers.extend(list(rev_primers.values()))
    elif isinstance(rev_primers, (list, tuple)):
        all_primers.extend(rev_primers)

    for primer in all_primers:
        binding = primer.binding.upper()
        if len(binding) < 10:
            continue
        rc_binding = get_rc(binding)

        for search_seq, strand in [(binding, 1), (rc_binding, -1)]:
            idx = 0
            while True:
                pos = seq_up.find(search_seq, idx)
                if pos == -1:
                    break
                key = (pos, pos + len(search_seq), strand)
                if key not in seen:
                    seen.add(key)
                    results.append({
                        "label": primer.name,
                        "start": pos,
                        "end":   pos + len(search_seq),
                        "type":  "Primer",
                        "strand": strand,
                        "full_sequence": primer.full_sequence,
                        "p_type": primer.p_type
                    })
                idx = pos + 1

    return results


def auto_annotate(sequence, global_features_dict):
    found_features = []
    seq_upper = sequence.upper()
    for key, feat in global_features_dict.items():
        feat_seq = feat.sequence.upper()
        if not feat_seq or len(feat_seq) < 5:
            continue
        start = 0
        while True:
            idx = seq_upper.find(feat_seq, start)
            if idx == -1:
                break
            found_features.append({"label": feat.label, "start": idx,
                                   "end": idx + len(feat_seq), "type": feat.type, "strand": 1})
            start = idx + 1
        feat_rc = get_rc(feat_seq)
        start = 0
        while True:
            idx = seq_upper.find(feat_rc, start)
            if idx == -1:
                break
            found_features.append({"label": f"{feat.label}(RC)", "start": idx,
                                   "end": idx + len(feat_rc), "type": feat.type, "strand": -1})
            start = idx + 1
    return found_features


def calculate_extension_time(bp_length, polymerase="Taq"):
    rate = POLYMERASE_RATES.get(polymerase, 60)
    return int((bp_length / 1000) * rate)


def adjust_features(features, offset, start_idx, end_idx, topology="Linear", template_len=0):
    res = []
    if not features:
        return []
    
    if topology == "Circular" and end_idx < start_idx:
        # Case where PCR product spans across the circular origin
        # Product is [start_idx : template_len] + [0 : end_idx]
        for f in features:
            s, e = f['start'], f['end']
            # Feature is in the first part [start_idx : template_len]
            if s >= start_idx:
                nf = f.copy()
                nf['start'] = s - start_idx + offset
                nf['end'] = e - start_idx + offset
                res.append(nf)
            # Feature is in the second part [0 : end_idx]
            elif e <= end_idx:
                nf = f.copy()
                nf['start'] = s + (template_len - start_idx) + offset
                nf['end'] = e + (template_len - start_idx) + offset
                res.append(nf)
    else:
        # Standard linear case
        for f in features:
            if f['start'] >= start_idx and f['end'] <= end_idx:
                nf = f.copy()
                nf['start'] = f['start'] - start_idx + offset
                nf['end'] = f['end'] - start_idx + offset
                res.append(nf)
    return res


def _hamming_find(text, pattern, max_mm=2):
    """Find first occurrence of pattern in text with at most max_mm mismatches.
    Returns start index or -1."""
    pl = len(pattern)
    if pl > len(text):
        return -1
    for i in range(len(text) - pl + 1):
        mm = sum(a != b for a, b in zip(text[i:i + pl], pattern))
        if mm <= max_mm:
            return i
    return -1


def _hamming_rfind(text, pattern, max_mm=2):
    """Find last occurrence of pattern in text with at most max_mm mismatches.
    Returns start index or -1."""
    pl = len(pattern)
    if pl > len(text):
        return -1
    result = -1
    for i in range(len(text) - pl + 1):
        mm = sum(a != b for a, b in zip(text[i:i + pl], pattern))
        if mm <= max_mm:
            result = i
    return result


def simulate_hr_logic_with_features(target_obj, amplicon_obj,
                                     max_end_skip=60, max_arm=600, min_arm=15,
                                     max_mismatches=2):
    """Insert를 target에 HR로 통합.
    - 양쪽 말단에서 max_end_skip bp까지는 비-homology 영역 허용.
    - max_mismatches: homology arm에서 허용하는 최대 미스매치 수 (기본 2bp).
    - homology arm 길이는 max_arm 부터 min_arm 까지 시도하며 가장 긴 매칭 우선.
    """
    t_seq, a_seq = target_obj.sequence, amplicon_obj.sequence
    a_len = len(a_seq)
    al_max = min(max_arm, a_len // 2)
    if al_max < min_arm:
        return None, [], 0, 0

    # --- Left arm 검색: offset 0..max_end_skip ---
    li, al, lo = -1, 0, 0
    end_skip = min(max_end_skip, a_len // 4)
    for offset in range(end_skip + 1):
        for l in range(al_max, min_arm - 1, -1):
            if offset + l > a_len:
                continue
            arm = a_seq[offset:offset + l]
            idx = t_seq.find(arm)
            if idx == -1 and max_mismatches > 0:
                idx = _hamming_find(t_seq, arm, max_mismatches)
            if idx != -1:
                li, al, lo = idx, l, offset
                break
        if li != -1:
            break
    if li == -1:
        return None, [], 0, 0

    # --- Right arm 검색: 끝에서 offset만큼 안쪽으로 들어가서 시도 ---
    ri, ar, ro = -1, 0, 0
    for offset in range(end_skip + 1):
        for l in range(al_max, min_arm - 1, -1):
            if offset + l > a_len:
                continue
            start_in_a = a_len - offset - l
            if start_in_a < lo + al:
                continue
            arm = a_seq[start_in_a:a_len - offset]
            idx = t_seq.rfind(arm)
            if idx == -1 and max_mismatches > 0:
                idx = _hamming_rfind(t_seq, arm, max_mismatches)
            if idx != -1 and idx >= li + al:
                ri, ar, ro = idx, l, offset
                break
        if ri != -1:
            break
    if ri == -1:
        return None, [], 0, 0

    # --- 통합 ---
    insert_part = a_seq[lo:a_len - ro]
    final_seq = t_seq[:li] + insert_part + t_seq[ri + ar:]

    # target 좌측 보존
    new_features = [f for f in target_obj.features if f['end'] <= li]
    # amplicon feature: insert_part 좌표계로 옮기고 target 위치(li)만큼 shift
    for f in amplicon_obj.features:
        fs, fe = f['start'], f['end']
        if fs >= lo and fe <= a_len - ro:
            nf = f.copy()
            nf['start'] = fs - lo + li
            nf['end'] = fe - lo + li
            new_features.append(nf)
    # target 우측 보존 (위치 shift)
    shift = len(insert_part) - (ri + ar - li)
    for f in target_obj.features:
        if f['start'] >= ri + ar:
            nf = f.copy()
            nf['start'] += shift
            nf['end'] += shift
            new_features.append(nf)

    return final_seq, new_features, li, li + len(insert_part)


# --- GFF/GBK parser ---

def _parse_gff_attributes(attr_str):
    out = {}
    for kv in attr_str.strip().split(';'):
        if '=' not in kv:
            continue
        k, v = kv.split('=', 1)
        out[k.strip()] = urllib.parse.unquote(v.strip())
    return out


def parse_gff3(file_path,
               wanted_types=("gene", "CDS", "tRNA", "rRNA", "ncRNA",
                             "snoRNA", "snRNA", "telomere", "centromere",
                             "ARS", "LTR_retrotransposon",
                             "transposable_element_gene")):
    features = []
    chromosomes = {}
    in_fasta = False
    current_chr = None
    current_chunks = []

    with open(file_path, 'r', encoding='utf-8', errors='replace') as f:
        for line in f:
            if not line:
                continue
            if line.startswith('##FASTA'):
                in_fasta = True
                continue
            if in_fasta:
                if line.startswith('>'):
                    if current_chr is not None:
                        chromosomes[current_chr] = ''.join(current_chunks)
                    current_chr = line[1:].strip().split()[0]
                    current_chunks = []
                else:
                    current_chunks.append(line.strip())
                continue
            if line.startswith('#') or not line.strip():
                continue
            parts = line.rstrip('\n').split('\t')
            if len(parts) < 9:
                continue
            seqid, source, ftype, start, end, score, strand, phase, attrs = parts[:9]
            if ftype not in wanted_types:
                continue
            try:
                start_i = int(start) - 1
                end_i = int(end)
            except ValueError:
                continue
            attr_dict = _parse_gff_attributes(attrs)
            label = attr_dict.get('gene') or attr_dict.get('Name') or attr_dict.get('ID') or ftype
            aliases = []
            for key in ('gene', 'Name', 'ID', 'Alias'):
                v = attr_dict.get(key, '')
                if v:
                    aliases.extend([x.strip() for x in v.split(',') if x.strip()])
            features.append({
                "label": label,
                "start": start_i,
                "end": end_i,
                "type": ftype,
                "strand": -1 if strand == '-' else 1,
                "chr": seqid,
                "attrs": attr_dict,
                "aliases": list(dict.fromkeys(aliases)),
            })
        if current_chr is not None:
            chromosomes[current_chr] = ''.join(current_chunks)

    if chromosomes:
        for f in features:
            chr_seq = chromosomes.get(f["chr"])
            if not chr_seq:
                continue
            sub = chr_seq[f["start"]:f["end"]]
            if f["strand"] == -1:
                sub = get_rc(sub)
            f["sequence"] = sub.upper()

    return features, chromosomes


def parse_annotation_file(file_path):
    features = []
    try:
        ext = os.path.splitext(file_path)[1].lower()
        if ext in ['.gbk', '.gb']:
            for record in SeqIO.parse(file_path, "genbank"):
                for feat in record.features:
                    label = feat.qualifiers.get('gene',
                                                 feat.qualifiers.get('locus_tag', [feat.type]))[0]
                    features.append({
                        "label": label,
                        "start": int(feat.location.start),
                        "end": int(feat.location.end),
                        "type": feat.type,
                        "strand": 1 if feat.location.strand in (1, None) else -1,
                        "chr": record.id,
                        "sequence": str(feat.extract(record.seq)),
                    })
        elif ext in ['.gff', '.gff3']:
            feats, _chrs = parse_gff3(file_path)
            features = feats
    except Exception:
        pass
    return features


def get_primer_analysis(seq, seq2=None):
    """
    primer3-py를 사용하여 프라이머의 물리적 성질 및 2차 구조 분석
    """
    import primer3
    from Bio.SeqUtils import molecular_weight, gc_fraction
    
    def analyze_single(s):
        s_upper = s.upper()
        mw = round(molecular_weight(s_upper, seq_type="DNA"), 2)
        gc = round(gc_fraction(s_upper) * 100, 1)
        tm = round(primer3.calc_tm(s_upper), 2)
        
        # Extinction Coefficient
        ext_coef = (s_upper.count('A')*15.4 + s_upper.count('C')*7.3 + s_upper.count('G')*11.7 + s_upper.count('T')*8.8) * 1000
        nmol_od = round(1000000 / ext_coef, 2) if ext_coef else 0
        ug_od = round(mw / (ext_coef/1000), 2) if ext_coef else 0
        
        hp = primer3.calc_hairpin(s_upper)
        hd = primer3.calc_homodimer(s_upper)
        
        return {
            "basic": {
                "Sequence": s_upper, "Length": len(s_upper), "MW": mw, 
                "GC%": gc, "Tm": tm, "nmole/OD": nmol_od, "ug/OD": ug_od
            },
            "hairpin": {
                "dg": round(hp.dg / 1000, 2), "tm": round(hp.tm, 2), 
                "dh": round(hp.dh / 1000, 2), "ds": round(hp.ds / 1000, 4),
                "warning": (hp.dg / 1000) < -3.0
            },
            "homodimer": {
                "dg": round(hd.dg / 1000, 2), 
                "dh": round(hd.dh / 1000, 2), "ds": round(hd.ds / 1000, 4),
                "structure": hd.ascii_structure if hasattr(hd, 'ascii_structure') else "N/A"
            }
        }

    results = {"p1": analyze_single(seq)}
    if seq2:
        results["p2"] = analyze_single(seq2)
        htd = primer3.calc_heterodimer(seq.upper(), seq2.upper())
        results["heterodimer"] = {
            "dg": round(htd.dg / 1000, 2),
            "dh": round(htd.dh / 1000, 2), "ds": round(htd.ds / 1000, 4),
            "structure": htd.ascii_structure if hasattr(htd, 'ascii_structure') else "N/A"
        }
    return results
