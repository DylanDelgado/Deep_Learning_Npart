#!/usr/bin/python

import random
import sys

random.seed()
num = int(sys.argv[1])
#print(random.randint(0, 2**29))
numbers = set();
for i in xrange(num):
    while True:
        x = random.randint(0, 2**29)
        if x not in numbers:
            numbers.add(x)
            break

#print((numbers))
output='\n'.join([str(i) for i in numbers])
#print(output)

fout = open('seedlist', 'w')
fout.write(output)
fout.write('\n')
fout.close()
#for i in numbers:
#	write(str(i))
#   write('\n')
