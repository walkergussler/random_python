def graphMeWithPlt(TransNet,fileList,output): #plots nx graph
    sparse_graph=csgraph_from_dense(TransNet)
    graph_obj=nx.from_scipy_sparse_matrix(sparse_graph,create_using=nx.DiGraph())
    nodeLabels={}
    edgeLabels={}
    for node in graph_obj.node:
    #------------------------------------------------------------------------------------------
    #this part generates names using some regex grab from your file names for cleaner charts
    #for no regex, use 
        # nodeLabels[node]=fileList[node]
        #or match the unique part of your filename and use that instead of the full name (below)
        # name=re.findall('.*_clipped/(.*)_unique.*',fileList[node])
        name=re.findall('(.*)_unique.*',fileList[node])
        # name=re.findall('(.*).*',fileList[node])
        nodeLabels[node]=name
    #------------------------------------------------------------------------------------------
    edgeListIterator=0
    for edge in sparse_graph.data:
        edgeLabels[graph_obj.edges()[edgeListIterator]]=str(edge)
        edgeListIterator+=1
    pos=nx.shell_layout(graph_obj,dim=2) 
    nx.draw_networkx_nodes(graph_obj,pos,node_shape="s")
    nx.draw_networkx_edges(graph_obj,pos,arrows=True)
    nx.draw_networkx_labels(graph_obj,pos,labels=nodeLabels)
    nx.draw_networkx_edge_labels(graph_obj,pos,edge_labels=edgeLabels)
    plt.axis('off')
    plt.savefig(output)

