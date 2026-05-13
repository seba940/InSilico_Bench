import json
import os
import traceback
from core.models import Primer, SequenceItem, Marker, EpitopeTag, Species, AnnotationRef, GlobalFeature

class LibraryManager:
    def __init__(self, filename="data/data.json"):
        self.filename = filename
        self.fwd_primers, self.rev_primers = {}, {}
        self.templates, self.amplicons, self.recombinants = {}, {}, {}
        self.digests   = {}  # restriction digest results
        self.ligations = {}  # ligation results
        self.markers, self.tags = {}, {}
        self.species, self.ann_refs = {}, {}
        self.global_features = {}
        self.settings = {}
        self.load()

    def save(self):
        data = {
            "species": [v.__dict__ for v in self.species.values()],
            "ann_refs": [v.__dict__ for v in self.ann_refs.values()],
            "primers": [v.__dict__ for v in list(self.fwd_primers.values()) + list(self.rev_primers.values())],
            "sequences": [v.__dict__ for v in
                          list(self.templates.values()) +
                          list(self.amplicons.values()) +
                          list(self.recombinants.values()) +
                          list(self.digests.values()) +
                          list(self.ligations.values())],
            "markers": [v.__dict__ for v in self.markers.values()],
            "tags": [v.__dict__ for v in self.tags.values()],
            "global_features": [v.__dict__ for v in self.global_features.values()],
            "settings": self.settings
        }
        with open(self.filename, 'w', encoding='utf-8') as f:
            json.dump(data, f, indent=4, ensure_ascii=False)

    def load(self):
        if not os.path.exists(self.filename):
            return
        try:
            with open(self.filename, 'r', encoding='utf-8') as f:
                d = json.load(f)
            self.settings = d.get("settings", {})
            for s in d.get("species", []):
                self.species[s['name']] = Species(**s)
            for a in d.get("ann_refs", []):
                self.ann_refs[a['name']] = AnnotationRef(**a)

            for p_data in d.get("primers", []):
                pt = p_data.get('p_type', p_data.get('type', 'Fwd'))
                p_data['p_type'] = pt
                if 'type' in p_data:
                    del p_data['type']
                p_obj = Primer(**{k: v for k, v in p_data.items()
                                  if k in Primer.__dataclass_fields__})
                if pt == "Fwd":
                    self.fwd_primers[p_obj.name] = p_obj
                else:
                    self.rev_primers[p_obj.name] = p_obj

            for s_data in d.get("sequences", []):
                cat = s_data.get('category', 'template')
                filtered = {k: v for k, v in s_data.items()
                            if k in SequenceItem.__dataclass_fields__}
                s_obj = SequenceItem(**filtered)
                if cat == "template" and not s_obj.kind:
                    s_obj.kind = "Genome" if s_obj.topology == "Linear" else "Plasmid"
                if cat == "template":
                    self.templates[s_obj.name] = s_obj
                elif cat == "amplicon":
                    self.amplicons[s_obj.name] = s_obj
                elif cat == "recombinant":
                    self.recombinants[s_obj.name] = s_obj
                elif cat == "digest":
                    self.digests[s_obj.name] = s_obj
                elif cat == "ligation":
                    self.ligations[s_obj.name] = s_obj

            self.global_features = {}

            for m_data in d.get("markers", []):
                self.markers[m_data['name']] = Marker(**m_data)
            for t_data in d.get("tags", []):
                self.tags[t_data['name']] = EpitopeTag(**t_data)

            def sync_internal_features(items_dict, prefix):
                for name, obj in items_dict.items():
                    main_type = "Marker" if prefix == "M" else "Tag"
                    key_full = f"{prefix}_{name}_Full"
                    self.global_features[key_full] = GlobalFeature(name, obj.sequence, main_type)
                    for fx in obj.features:
                        label = fx['label']
                        seq_sub = obj.sequence[fx['start']:fx['end']]
                        if len(seq_sub) > 5:
                            key_int = f"Internal_{label}"
                            self.global_features[key_int] = GlobalFeature(label, seq_sub, fx['type'])

            sync_internal_features(self.markers, "M")
            sync_internal_features(self.tags, "T")

            for g in d.get("global_features", []):
                if g['label'] not in self.global_features:
                    self.global_features[g['label']] = GlobalFeature(**g)

        except Exception:
            print(f"Load Error: {traceback.format_exc()}")
