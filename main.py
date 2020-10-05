import table_transformer as tt
import consensus
import time
import branch_list
start_time = time.time()

df = tt.TableTransformer('/home/masha/mut_exp.tsv', "/home/masha/tree/Listeria_monocytogenes_68_83/mutations/consensuses.fasta")

list_bunch = branch_list.MakeBranchList('/home/masha/tree/trees/nj-global-tree.tre')[0]
single_bunch_set = branch_list.MakeBranchList('/home/masha/tree/trees/nj-global-tree.tre')[1]

if __name__ == '__main__':
    #print(df.group_mutation(list_bunch, single_bunch_set))
    df.group_mutation(list_bunch, single_bunch_set).to_csv('MutationSignature/out.csv')

print("--- %s seconds ---" % (time.time() - start_time))
# See PyCharm help at https://www.jetbrains.com/help/pycharm/
