from Bio import SeqIO

class Consensus(object):
    def __init__(self, file_path):
        self.block_sequence = SeqIO.to_dict(SeqIO.parse(file_path, "fasta"))

    def get_consensus(self, block_names, position):
        return self.block_sequence[block_names].seq[position]
    #def get_len_block(self, block_names):
    #    return len(self.block_sequence[block_names])
