import webbrowser

class BaseProvider:
    name = "Generic"
    base_url = ""
    
    def get_locus_url(self, ident):
        return None

    def get_sequence_url(self, ident):
        return None

class SGDProvider(BaseProvider):
    name = "SGD"
    base_url = "https://www.yeastgenome.org/locus/"
    
    def get_locus_url(self, ident):
        return self.base_url + ident

    def get_sequence_url(self, ident):
        return self.base_url + ident + "#sequence"

class NCBIProvider(BaseProvider):
    name = "NCBI"
    base_url = "https://www.ncbi.nlm.nih.gov/gene/?term="
    
    def get_locus_url(self, ident):
        return self.base_url + ident

    def get_sequence_url(self, ident):
        return f"https://www.ncbi.nlm.nih.gov/nucleotide/?term={ident}"

class PomBaseProvider(BaseProvider):
    name = "PomBase"
    base_url = "https://www.pombase.org/gene/"
    
    def get_locus_url(self, ident):
        return self.base_url + ident

    def get_sequence_url(self, ident):
        return self.base_url + ident

class FlyBaseProvider(BaseProvider):
    name = "FlyBase"
    base_url = "https://flybase.org/reports/"
    
    def get_locus_url(self, ident):
        return self.base_url + ident

class WormBaseProvider(BaseProvider):
    name = "WormBase"
    base_url = "https://wormbase.org/species/all/gene/"
    
    def get_locus_url(self, ident):
        return self.base_url + ident

class ProviderRegistry:
    _providers = {
        "Saccharomyces cerevisiae": SGDProvider(),
        "Homo sapiens": NCBIProvider(),
        "Mus musculus": NCBIProvider(),
        "Schizosaccharomyces pombe": PomBaseProvider(),
        "Drosophila melanogaster": FlyBaseProvider(),
        "Caenorhabditis elegans": WormBaseProvider(),
        "Generic": NCBIProvider()
    }
    
    _db_name_map = {
        "SGD": SGDProvider(),
        "NCBI": NCBIProvider(),
        "PomBase": PomBaseProvider(),
        "FlyBase": FlyBaseProvider(),
        "WormBase": WormBaseProvider()
    }

    @classmethod
    def get_by_species(cls, species_name):
        return cls._providers.get(species_name, cls._providers["Generic"])

    @classmethod
    def get_by_db_name(cls, db_name):
        return cls._db_name_map.get(db_name, cls._db_name_map["NCBI"])

    @classmethod
    def list_species(cls):
        return sorted(list(cls._providers.keys()))
