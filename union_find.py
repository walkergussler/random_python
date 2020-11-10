from __future__ import absolute_import
from __future__ import division
from __future__ import print_function

import numpy as np
from scipy.sparse import dok_matrix

class UnionFind(object):
  # An implementation of union find data structure.
  # It uses weighted quick union by rank with path compression.
  def __init__(self, node_ids):
    # Initialize an empty union find object with N items.
    #  N: Number of items in the union find object.
    self._sets = {
        node_id : {
          'rank' : 0,
          'parent' : node_id
        } for idx,node_id in enumerate(node_ids)
    }

  def find(self, x):
    try:
      p_idx = self._sets[x]['parent']
      if p_idx != x:
        self._sets[x]['parent'] = self.find(self._sets[p_idx]['parent'])

      return self._sets[x]['parent']

    except KeyError:
      raise KeyError('ID {0} is not a member of the union'.format(x))

  def join(self, p, q):
    # Combine sets containing p and q into a single set.
    p_id = self.find(p)
    q_id = self.find(q)
    pRoot = self._sets[p_id]
    qRoot = self._sets[q_id]

    if pRoot['rank'] < qRoot['rank']:
      pRoot['parent'] = q_id

    elif pRoot['rank'] > qRoot['rank']:
      qRoot['parent'] = p_id

    else:
      qRoot['parent'] = p_id
      pRoot['rank'] += 1

  def connected(self):
    it = iter(self._sets)
    f = self.find(next(it))
    for i in it:
      if self.find(i) != f:
        return False

    return True
