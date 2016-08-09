import sys

from py_rstr_max.rstr_max import *

def produce_seeds(strings, min_size=1, max_size=9000):
  rstr = Rstr_max()
  for string in strings:
    rstr.add_str(string)
  r = rstr.go()

  stack = []
  for (offset_end, nb), (l, start_plage) in r.iteritems():
    seed = rstr.global_suffix[offset_end-l:offset_end]
    if len(seed) < min_size or len(seed) >  max_size or ("N" in seed):
      continue
    stack.append(tuple((nb,seed)))
  stack.sort(key=lambda x:x[0], reverse=True)
  return stack
