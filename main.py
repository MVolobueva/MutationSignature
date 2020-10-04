import table_transformer as tt
import consensus
import time
start_time = time.time()

df = tt.TableTransformer('/home/masha/tree/Listeria_monocytogenes_68_83/mutations/mut.tsv', "/home/masha/tree/Listeria_monocytogenes_68_83/mutations/consensuses.fasta")
with open('/home/masha/tree/trees/nj-global-tree.tre',  'r') as f:
    a = 0
    line = [i.strip() for i in f.readlines() ]
    for i in range(len(line)):
        if ')' in line[i]:
            line[i] = ')'
    bunch = set()
    list_bunch = []
    single_bunch_set = set()
    while ')' in line:
        d = (line.index(')'))
        line.pop(d)
        while line[d-1] != '(':
            single_bunch_set.add(line[d-1].split(':')[0])
            bunch.add(line[d-1].split(':')[0])
            d -= 1
        line.pop(d-1)
        a += 1
        list_bunch.append(bunch)
        bunch = set()
    ls = [{i} for i in single_bunch_set]
    list_bunch = list_bunch + ls
print(list_bunch)
# Press the green button in the gutter to run the script.
if __name__ == '__main__':
    #print(df.group_mutation(list_bunch, single_bunch_set))
    df.group_mutation(list_bunch, single_bunch_set).to_csv('/home/masha/PycharmProjects/mutation_signature/out.csv')

print("--- %s seconds ---" % (time.time() - start_time))
# See PyCharm help at https://www.jetbrains.com/help/pycharm/
