import pandas as pd
import numpy as np
import consensus
from Bio import SeqIO

class TableTransformer(object):
    def __init__(self, file_path, path_to_consensus):
        self.table = pd.read_csv(file_path, sep='\t')
        self.block_sequence = SeqIO.to_dict(SeqIO.parse(path_to_consensus, "fasta"))

    def _get_consensus(self, block_names, position):
        return self.block_sequence[block_names].seq[position]

    def _drop_none(self):
        self.table = self.table.replace('.', np.nan)
        self.table = self.table.fillna(method='pad')
        return self.table

    def _add_consensus(self):
        self.table = self.table[self.table['stop(gaps)/change'].str.isalpha()]
        self.table['consensus'] = self.table.apply(lambda x: self._get_consensus(x['block'], x['start']), axis = 1)

        return self.table

    def _organism_name(self):
        self.table['organism_name'] = self.table.apply(lambda x: x['fragment'].split('&')[0], axis = 1)
        return self.table

    def group_mutation(self, clade_list, single_bunch_set):
        self._drop_none()
        self._add_consensus()
        self._organism_name()
        df1 = self.table.set_index(['block', 'consensus', 'stop(gaps)/change', 'start'])
        #print(list(set(df1.index))[1])
        index_ls = set(df1.index)
        #print([set(df1.loc[i]['organism_name']) for i in index_ls])

        self.table = pd.DataFrame({'block': [i[0] for i in index_ls], \
                                    'consensus': [i[1] for i in index_ls], \
                                    'change': [i[2] for i in index_ls],\
                                    'start': [i[3] for i in index_ls],\
                                    'set of clade': [set(df1.loc[i]['organism_name']) for i in index_ls]})
        self.table['is clade'] = self.table['set of clade'].apply(lambda x: x in clade_list or x in single_bunch_set - x)
        return self.table

