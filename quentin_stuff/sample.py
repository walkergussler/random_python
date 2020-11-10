class quasiMedianAlignment(Mapping):
    def __init__(self, ali):
        self.ali = ali
        self.qm = {}
        self.qm_reverse_index = {}

    def __iter__(self):
        for i in self.ali:
            yield i

        for i in self.qm:
            yield i

    def __len__(self):
        return len(self.ali) + len(self.qm)

    def __getitem__(self, key):
        try:
            return self.ali.get(key)
        except KeyError:
            return self.qm[key]

    def __delitem__(self, key):
        del self.qm[key]
        del self.qm_reverse_index[key]

    def hamming(self, p, q):
        p_i = self[p]
        q_i = self[q]
        diff = p_i != q_i
        return self.ali.col_counts[diff].sum()

def spanningNetwork(ali,g):
    dist = [ (ali.hamming(i,j),(i,j)) for j,i in combinations(ali, 2) ]
    numNodes=len(ali)
    gval=int(math.ceil(85*(-math.exp(-numNodes/1000)+1)))
    gstr='gray'+str(gval)
    # G = nx.Graph()
    # G.add_nodes_from(ali)
    uf = UnionFind(ali)
    heapq.heapify(dist) 
    while not uf.connected():
        current,_ = dist[0]
        while dist[0][0] == current:
            d,(i,j) = heapq.heappop(dist)
            if d>=10:
                g.add_edge(i,j,len=10,color="red",label=d,fontsize="24")
            else:
                g.add_edge(i,j,len=d,color=gstr)
            uf.join(i,j)      
    # max_weight = current
    # d = max_weight
    # while d < max_weight + delta:
    # d,(i,j) = heapq.heappop(dist)
    # G.add_edge(i,j, dist=d)
    return g
