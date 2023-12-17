def load_queries(path):
    with open(path,"r") as f:
        queries=[line.rstrip() for line in f.readlines()]
    return queries


if __name__=="__main__":
    queries=load_queries("query.fasta")
    print(queries[0])