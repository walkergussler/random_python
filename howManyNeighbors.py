from __future__ import division
import math
from decimal import *
import sys

def howManyNeighbors(radius,seqlen=264,alphabet=4):
    neighbors=0
    for a in range(radius):
        x=a+1
        constant=Decimal((alphabet-1)**x/math.factorial(x))
        variable=1
        for i in range(x):
            term=seqlen-i
            variable=Decimal(term*variable)
        neighbors+=Decimal(constant*variable)
    return neighbors
    
# print(howManyNeighbors(20,264,4))
# print(howManyNeighbors(40,264,4))
# print(howManyNeighbors(60,264,4))
# print(howManyNeighbors(80,264,4))
# print(howManyNeighbors(100,264,4))
# print(howManyNeighbors(120,264,4))
# print(howManyNeighbors(264,264,4))
# print('%e' % 4**264)
# print("radius, seqlen, alphabet")
# print(howManyNeighbors(3))
# print(howManyNeighbors(3,264,4))
print(howManyNeighbors(int(sys.argv[1])))