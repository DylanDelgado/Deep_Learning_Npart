#!/usr/bin/env python3

import random
import sys

random.seed()
num = int(sys.argv[1])

no_fp = False
try:
    fp = open('rdlist')
except IOError:
    no_fp = True

oldnum = []
if not no_fp:
    with fp:
        for line in fp:
            oldnum.append(int(line))

numbers = []
seen = set(oldnum)
for i in range(num):
    while True:
        x = random.randint(0, 2**29)
        if x not in seen:
            numbers.append(x)
            seen.add(x)
            break

output='\n'.join([str(i) for i in numbers])
print(output)
