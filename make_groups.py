        else:
            h.add_node(f1)
            h.add_node(f2)
    cc=0
    gdic={}
    hdic={}
    for comp in nx.connected_components(g):
        cc+=1
        for file in comp:
            gdic[file]=str(cc)
    ccc=0
    for comp in nx.connected_components(h):
        ccc+=1
        for file in comp:
            gdic[file]=str(ccc)
    
    for file in set(files):
        
    for file in set(files):
        ','.join([file,gdic[file],hdic[file])
