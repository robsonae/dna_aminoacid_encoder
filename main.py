
codons = {
    'ATA': 'I', 'ATC': 'I', 'ATT': 'I', 'ATG': 'M', 'ACA': 'T',
    'ACC': 'T', 'ACG': 'T', 'ACT': 'T', 'AAC': 'N', 'AAT': 'N',
    'AAA': 'K', 'AAG': 'K', 'AGC': 'S', 'AGT': 'S', 'AGA': 'R',
    'AGG': 'R', 'CTA': 'L', 'CTC': 'L', 'CTG': 'L', 'CTT': 'L',
    'CCA': 'P', 'CCC': 'P', 'CCG': 'P', 'CCT': 'P', 'CAC': 'H',
    'CAT': 'H', 'CAA': 'Q', 'CAG': 'Q', 'CGA': 'R', 'CGC': 'R',
    'CGG': 'R', 'CGT': 'R', 'GTA': 'V', 'GTC': 'V', 'GTG': 'V',
    'GTT': 'V', 'GCA': 'A', 'GCC': 'A', 'GCG': 'A', 'GCT': 'A',
    'GAC': 'D', 'GAT': 'D', 'GAA': 'E', 'GAG': 'E', 'GGA': 'G',
    'GGC': 'G', 'GGG': 'G', 'GGT': 'G', 'TCA': 'S', 'TCC': 'S',
    'TCG': 'S', 'TCT': 'S', 'TTC': 'F', 'TTT': 'F', 'TTA': 'L',
    'TTG': 'L', 'TAC': 'Y', 'TAT': 'Y', 'TAA': '*', 'TAG': '*',
    'TGC': 'C', 'TGT': 'C', 'TGA': '*', 'TGG': 'W'
}

protein_weights = {
    'A': 89.0932,
    'C': 121.1582,
    'D': 133.1027,
    'E': 147.1293,
    'F': 165.1891,
    'G': 75.0666,
    'H': 155.1546,
    'I': 131.1729,
    'K': 146.1876,
    'L': 131.1729,
    'M': 149.2113,
    'N': 132.1179,
    'O': 255.3134,
    'P': 115.1305,
    'Q': 146.1445,
    'R': 174.201,
    'S': 105.0926,
    'T': 119.1192,
    'U': 168.0532,
    'V': 117.1463,
    'W': 204.2252,
    'Y': 181.1885,
    '*': 0.0,
    'X': 0.0
}

def reverse_complement(dna: str) -> str:
    comp = str.maketrans("ACGTNacgtn", "TGCANtgcan")        # N -> niejednoznaczne kodony, błędy
    return dna.translate(comp)[::-1]

def translate_frame(seq: str, offset: int) -> str:
    pep = []
    for i in range(offset, len(seq) -2, 3):
        codon = seq[i:i+3]
        pep.append(codons.get(codon, 'X'))          #accesses codons dict, if the key is not there the value 'X' is returned
    return ''.join(pep)

class Sequence:

    def __init__(self, sequence: str):

        self.sequence_raw = ''.join(sequence.split()).upper()                       # surowa sekwencja - wszystkie symbole
        self.sequence = ''.join(ch for ch in self.sequence_raw if ch in 'ATGC')     # wyczyszczona sekwencja - tylko ATGC

        n = len(sequence) or 1          # unik podziału przez 0
        self.a_percentage = self.sequence.count('A') / n * 100
        self.t_percentage = self.sequence.count('T') / n * 100
        self.g_percentage = self.sequence.count('G') / n * 100
        self.c_percentage = self.sequence.count('C') / n * 100

        # 3 forward open reading frames
        self.orf1transl = translate_frame(self.sequence, 0)
        self.orf2transl = translate_frame(self.sequence, 1)
        self.orf3transl = translate_frame(self.sequence, 2)

        # 3 reverse open reading frames (nić komplementarna)
        rc = reverse_complement(self.sequence)
        self.orf4transl = translate_frame(rc, 0)
        self.orf5transl = translate_frame(rc, 1)
        self.orf6transl = translate_frame(rc, 2)

        #przeliczenie mas białek
        self.prot_weight_orf1 = self.prot_seq_weight(self.orf1transl)
        self.prot_weight_orf2 = self.prot_seq_weight(self.orf2transl)
        self.prot_weight_orf3 = self.prot_seq_weight(self.orf3transl)
        self.prot_weight_orf4 = self.prot_seq_weight(self.orf4transl)
        self.prot_weight_orf5 = self.prot_seq_weight(self.orf5transl)
        self.prot_weight_orf6 = self.prot_seq_weight(self.orf6transl)

    def prot_seq_weight(self, prot_seq: str) -> float:
        return sum(protein_weights.get(aa, 0.0) for aa in prot_seq if aa != '*')

def read_fasta(path: str) -> str:
    seq_parts = []
    with open(path, 'r', encoding="utf-8") as fin:
        for line in fin:
            if not line:        # an empty line has a False value in Python
                continue
            if line.startswith('>'):
                if seq_parts:       # if > is encountered and a sequence has already been added to seq_parts -> break...
                                    # ...only one fasta sequence per file is added. safeguard for more sequences in a file.
                    break
                else:
                    continue
            seq_parts.append(line.strip())
    return ''.join(seq_parts)

from textwrap import fill as _fill

def wrap120(text: str) -> str:          # code to wrap text to output file at 60 chars
    return _fill(text, width=120)

if __name__== "__main__":

    in_path = 'sequence_in.fasta'
    out_path = 'sequence_out.txt'

    seq_str = read_fasta(in_path)
    s = Sequence(seq_str)

    ncl_content_info = (
        "Nucleotide content (only A/C/G/T counted):\n"
        f" Adenine: {s.a_percentage:.2f}%\t"
        f" Thymine: {s.t_percentage:.2f}%\t"
        f" Guanine: {s.g_percentage:.2f}%\t"
        f" Cytosine: {s.c_percentage:.2f}%\n\n"
    )

    def info_block(title: str, pep: str, mass: float) -> str:
        return (
            f"{title}\n"
            f"{wrap120(pep)}\n"
            f"Molecular weight: {mass:.2f} g/mol\n\n"
        )

    infoORF1 = info_block("ORF1 (forward) amino acid sequence: ", s.orf1transl, s.prot_weight_orf1)
    infoORF2 = info_block("ORF2 (forward) amino acid sequence ", s.orf2transl, s.prot_weight_orf2)
    infoORF3 = info_block("ORF3 (forward) amino acid sequence:", s.orf3transl, s.prot_weight_orf3)
    infoORF4 = info_block("ORF4 (reverse) amino acid sequence:", s.orf4transl, s.prot_weight_orf4)
    infoORF5 = info_block("ORF5 (reverse) amino acid sequence:", s.orf5transl, s.prot_weight_orf5)
    infoORF6 = info_block("ORF6 (reverse) amino acid sequence:", s.orf6transl, s.prot_weight_orf6)

    with open(out_path, 'w', encoding='utf-8') as fout:
        fout.write(ncl_content_info)
        fout.write(infoORF1)
        fout.write(infoORF2)
        fout.write(infoORF3)
        fout.write(infoORF4)
        fout.write(infoORF5)
        fout.write(infoORF6)

    print(f"Wyniki konwersji zapisano do pliku: {out_path}")