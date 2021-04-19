
from bio_structures import DNA_Codons, Nucleotides
import random
from collections import Counter


class bio_seq():

    def __init__(self, seq="ATCG", seq_type='DNA', label = 'No label'):
        self.seq = seq.upper()
        self.label = label
        self.seq_type = seq_type
        self.is_valid = self.__validate()
        assert self.is_valid, f"provided data doesn't seem to be correct {self.seq_type}sequence"

    def __validate(self):
        #Private function; does not appear in the drop manu
       return set(Nucleotides).issuperset(self.seq)

    def get_seq_type(self):
        return f'Sequence type {self.seq_type}'

    def get_seq_info(self):
        """Returns 4 strings; full sequence information"""
        return f'[Label]: {self.label}\n [Sequence]: {self.seq}\n [Type]: {self.seq_type}\n [Length]: {len(self.seq)}'

    def gen_random_seq(self, length=10, seq_type='DNA'):
        """Generates a random DNA sequence"""
        seq = ''.join([random.choice(Nucleotides) for x in range(length)])
        self.__init__(seq, seq_type, 'Randomly generated sequence')
        # Reinitialize to check if the sequence is valid

    def count_nuc_freq(self):
        return dict(Counter(self.seq))

    def transcribe(self):
        return self.seq.replace("T", "U")

    def reverse_complement(self):
        # Swapping Adenine with Thymine and guanine with cytozine and vice versa.
        mapping = str.maketrans('ATCG', 'TAGC')
        return self.seq.translate(mapping)[::-1]

    def gc_content(self):
        return round((self.seq.count('G') + self.seq.count('C')) / len(self.seq) * 100)

    def gc_in_subsec(self, k=20):
        res = []
        for i in range(0, len(self.seq) - k + 1, k):
            sub_seq = self.seq[i:i + k]
            res.append(round((sub_seq.count('G') + sub_seq.count('C')) / len(sub_seq) * 100))
        return res

    def translate(self, init_pos=0):
        return [DNA_Codons[self.seq[pos:pos + 3]] for pos in range(init_pos, len(self.seq) - 2, 3)]

    def codon_usage(self, amino_acid):
        tmp_list = []
        for i in range(0, len(self.seq) - 2, 3):
            if DNA_Codons[self.seq[i:i + 3]] == amino_acid:
                tmp_list.append(self.seq[i:i + 3])

        freq_dict = dict(Counter(tmp_list))
        total = sum(freq_dict.values())
        for seq in freq_dict:
            freq_dict[seq] = round(freq_dict[seq] / total, 2)
        return freq_dict

    def gen_read_frames(self):
        tmp_seq = bio_seq(self.reverse_complement(), self.seq_type)
        frames = [self.translate(0), self.translate(1), self.translate(2),
                  tmp_seq.translate(0), tmp_seq.translate(1), tmp_seq.translate(2)]
        del tmp_seq
        return frames

    def proteins_from_rf(self, aa_seq):
        current_protein = []
        proteins = []
        for aa in aa_seq:
            if aa == "_":
                # Stops accumulating amino acid if STOP codon was found
                if current_protein:
                    for p in current_protein:
                        proteins.append(p)
                        current_protein = []
            else:
                # Starts accumulating aa if M start codon was found:
                if aa == "M":
                    current_protein.append("")
                for i in range(len(current_protein)):
                    current_protein[i] += aa
        return proteins

    def all_proteins_from_RF(self, start_pos=0, end_pos=0, ordered=False):
        # Computes all possible proteins from open reading frames
        if end_pos > start_pos:
            # this is an option to check a specific range in a seq with defined start and end position
            tmp_Seq = bio_seq(self.seq[start_pos: end_pos], self.seq_type)
            rfs = tmp_Seq.gen_read_frames()
        else:
            # or just check the whole sequence without defining start and end position
            rfs = self.gen_read_frames()

        res = []
        for rf in rfs:
            proteins = self.proteins_from_rf(rf)
            for p in proteins:
                res.append(p)

        if ordered:
            # Sorted option is added by length or it can be alphabetical
            return sorted(res, key=len, reverse=True)
            # reverse parameter returns the protein list from the longest protein
        return res