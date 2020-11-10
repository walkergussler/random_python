import csv
import itertools
print("--------------------")
with open('Indiana.csv') as csvfile:
    theList=[]
    eqcount=0
    r=csv.reader(csvfile)
    for line in r:
        if float(line[2])<.037:
            theList.append(line)
    theLists=[[theList[0][0],theList[0][1]],[theList[3][0],theList[3][1]]]
    biglist=[]
    for line in theList:  
        source=line[0]
        target=line[1]
        biglist.append(source)
        biglist.append(target)
        found=False
        for id in range(len(theLists)):
            if source in theLists[id] or target in theLists[id]:
                theLists[id].extend([source,target])
                found=True
        if found==False:
            theLists.append([source,target])
    theset=[]
    bigset=set(biglist)
    d={}
    for item in bigset:
        d[item]=0
    print(len(bigset))
    for i in theLists:
        if i not in theset:
            theset.append(set(i))
    finalList=[]
    for item in theset:
        # print(list(item))
        finalList.append(list(item))
    editedList=[]
    for item in finalList:
        print(item)
        # for subitem in item:
            # for key in d:
                # if key==subitem:
                    # d[key]+=1
                    # if d[key]==1:
                        # newList.append(key)
                    # else:
                        # print(d[key])
        # editedList.append(newList)
    # flatlist=[]
    # for k in editedList:
        # if k==[]:
            # editedList.remove(k)
    # for item in editedList:
        # print(item)
    