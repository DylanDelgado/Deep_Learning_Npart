#!/usr/bin/python

import random
import sys

random.seed()
num = int(sys.argv[1])
#print(random.randint(0, 2**29))
numbers = []
for i in xrange(num):
	while True:
		x = random.randint(0, 2**29)
		if x not in numbers:
			numbers.append(x)
			break

#output=' '.join([str(i) for i in numbers])
#print(output)

for i in numbers:
	print(i)
