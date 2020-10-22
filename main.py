import table_transformer as tt
import consensus
import time
import branch_list
import os
start_time = time.time()
for i in os.listdir('example'):
    df = tt.TableTransformer('example/' + i +'/mut.tsv', "example/" + i + "/consensuses.fasta")
    list_bunch = branch_list.MakeBranchList('example/' + i+ '/nj-global-tree.tre')[0]
    single_bunch_set = branch_list.MakeBranchList('example/' + i+ '/nj-global-tree.tre')[1]

    if __name__ == '__main__':
        #print(df.group_mutation(list_bunch, single_bunch_set))
        df.group_mutation(list_bunch, single_bunch_set)[1].to_csv('example/' + i +'/out.csv')

print("--- %s seconds ---" % (time.time() - start_time))
# See PyCharm help at https://www.jetbrains.com/help/pycharm/
