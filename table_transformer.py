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

    def _get_len_block(self, block_names):
        return len(self.block_sequence[block_names])

    def _drop_none(self):
        self.table = self.table.replace('.', np.nan)
        self.table = self.table.fillna(method='pad')
        return self.table

    def _add_consensus(self):
        self.table = self.table[self.table['stop(gaps)/change'].str.isalpha()]
        self.table['consensus'] = self.table.apply(lambda x: self._get_consensus(x['block'], x['start']), axis = 1)

        return self.table

    def _add_block_len(self):
        self.table['block_len'] = self.table.apply(lambda x: int(x['block'].split('x')[1].split('n')[0]), axis = 1)
        return self.table

    def _organism_name(self):
        self.table['organism_name'] = self.table.apply(lambda x: x['fragment'].split('&')[0], axis = 1)
        return self.table

    def _filter_by_block_name(self):
        self.table['block_name'] = self.table.apply(lambda x: x['block'][0], axis = 1)

        #self.table = self.table[self.table['block_name'] == 's']
        return self.table[self.table['block_name'] == 'r']

    def group_mutation(self, clade_list, single_bunch_set):
        self._drop_none()
        self._filter_by_block_name()
        self._add_block_len()
        self._add_consensus()
        self._organism_name()
        df1 = self.table.set_index(['block', 'consensus', 'stop(gaps)/change', 'start', 'block_len'])
        index_ls = set(df1.index)


        self.table = pd.DataFrame({'block': [i[0] for i in index_ls], \
                                    'consensus': [i[1] for i in index_ls], \
                                    'change': [i[2] for i in index_ls],\
                                    'start': [i[3] for i in index_ls],\
                                    'block_len': [i[4] for i in index_ls],\
                                    'set of clade': [set(df1.loc[i]['organism_name']) for i in index_ls]})
        self.table = self.table.sort_values(['block', 'start'], ascending=[True, True])
        self.table['start1'] = self.table['start'].diff()
        self.table['start2'] = self.table['start'].diff(-1)
        self.table = self.table[self.table['block_len'] - self.table['start'] > 2]
        self.table = self.table[self.table['start1'] != 2 ]
        self.table = self.table[self.table['start1'] != 1]
        self.table = self.table[self.table['start2'] != -2]
        self.table = self.table[self.table['start2'] != -1]


        self.table['is clade'] = self.table['set of clade'].apply(lambda x: x in clade_list or x in single_bunch_set - x)
        self.table = self.table[self.table['is clade']]
        dt = {'T>G':'A>C', "T>C":'A>G', 'T>A':'A>T', 'G>T':'C>A', 'G>C':'C>G', 'G>A':'C>T'}
        self.table['mutation'] = self.table.apply(lambda x: x['consensus'] + '>' + x['change'], axis = 1)

        self.table['mut'] = self.table.apply(lambda x: dt[x['mutation']] if x['mutation'] in dt.keys() else x['mutation'], axis = 1)
        return self.table, self.table['mut'].value_counts().sort_values()

