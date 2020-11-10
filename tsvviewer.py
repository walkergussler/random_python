#!/usr/bin/env python3
import sys
import os

infile=sys.argv[1]

with open(infile) as f:
    lines=f.readlines()
for line in lines:
    splitline=line.split("\t")
    fmt='{:>12} '*len(splitline)
    stripped=fmt.rstrip()
    cleansplit=
    print(splitline)
    print(stripped)
    newline=stripped.format(splitline)
    print(newline)

