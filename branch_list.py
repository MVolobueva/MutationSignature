def MakeBranchList(file_path):
    with open(file_path, 'r') as f:
        line = [i.strip() for i in f.readlines()]
        for i in range(len(line)):
            if ')' in line[i]:
                line[i] = ')'
        bunch = set()
        list_bunch = []
        single_bunch_set = set()
        while ')' in line:
            d = (line.index(')'))
            line.pop(d)
            while line[d - 1] != '(':
                single_bunch_set.add(line[d - 1].split(':')[0])
                bunch.add(line[d - 1].split(':')[0])
                d -= 1
            line.pop(d - 1)
            list_bunch.append(bunch)
            bunch = set()
        ls = [{i} for i in single_bunch_set]
        list_bunch = list_bunch + ls
    return list_bunch, single_bunch_set
