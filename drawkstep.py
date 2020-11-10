from __future__ import absolute_import
from __future__ import print_function
from __future__ import division

"""\
Draw a k-step network from sequences.

..moduleauthor: Joseph (Walker) Gussler <mnz0@cdc.gov>
..moduleauthor: Seth Sims <xzy3@cdc.gov>
"""

import argparse
import sys
import os.path
import math
import subprocess
import warnings
from itertools import ifilter,imap,count,chain
from collections import defaultdict,Counter
from tempfile import NamedTemporaryFile
from operator import itemgetter

import six
import networkx as nx
from Bio import SeqIO
from toolz.itertoolz import first

from ghost.util import msa
from ghost.util.distance import chunked_hamming

from filigree.network import kstep,contract_short_edges
from filigree.util.jinja2 import render_to_file
import filigree.util.logging
filigree.util.logging.consoleLogConfig()

DEFAULT_LEGEND_FONT_SIZE = 96
LEGEND_FONT_SIZE = {
  4 : 8,
  10 : 12,
  25 : 24,
  80 : 70
}

def determineLegendFontSize(g):
  lsp = max(l for l in chain.from_iterable(six.itervalues(nx.all_pairs_dijkstra_path_length(g, weight='len'))))
  fontSize=DEFAULT_LEGEND_FONT_SIZE
  for len_bin, font_size in sorted(LEGEND_FONT_SIZE.items(), key=lambda t: t[0]):
    if lsp < len_bin:
      fontSize = font_size
      break

  return str(fontSize)

def add_edge_attributes(g, drawmode):
  gstr = 'gray{}'.format(int(math.ceil(70*(-math.exp(-g.number_of_nodes()/1000)+1)+20)))
  for u, v, data in g.edges_iter(data=True):
    dist = data['len']
    if dist < 10:
      data['color'] = gstr
      if drawmode=='weights':
        data['label']=dist
    else:
      assert drawmode != 'dont', 'programming error kstep mode dont edge above threshold'

      if drawmode=='white':
        color1=g.node[u]['color']
        color2=g.node[v]['color']
        if color1!=color2 and color1!=0 and color2 !=0:
          data['color'] = 'white'
          data['len'] = 9
      elif drawmode=='red':
        color1=g.node[u]['color']
        color2=g.node[v]['color']
        if color1!=color2 and color1!=0 and color2 !=0:
          data['color'] = 'red'
          data['len'] = 9
          data['label'] = dist
      elif drawmode=='allred':
        data['color'] = 'red'
        data['len'] = 9
        data['label'] = dist
      elif drawmode=='weights':
        data['color'] = gstr
        data['label'] = dist
        data['len'] = 9
      else:
        raise ValueError("You have entered an invalid option for drawmode. Look at help for more info")

  return g

def draw_kstep(g, dist_iter, draw_mode, output_name, keep_tmp_files=False, output_type='png', color_scheme='set28'):
  threshold = float('inf')
  if draw_mode == 'dont':
    threshold = 10

  sample_names = set(data['sample'] for _,data in g.nodes_iter(data=True))
  kstep(dist_iter, g, threshold=threshold)

  def merge_node(g,u,v):
    if g.node[u]['sample'] != g.node[v]['sample']:
      sample_names.add('shared')
      g.node[u]['sample'] = 'shared'
      g.node[v]['sample'] = 'shared'

  contract_short_edges(g, 1, merge_node=merge_node)

  color_id = count(2)
  colors = defaultdict(lambda: str(next(color_id)))
  colors['shared'] = '1'

  numNodes = g.number_of_nodes()
  w_i = .12
  if numNodes <= 3000:
    w_i = math.exp(-numNodes/500)*.2+.1

  set_seq_count = Counter()
  for _,data in g.nodes_iter(data=True):
    freq = data['freq']
    s = data.get('sample', '_no_sample_id_')
    if 'color' not in data:
      data['color'] = colors[s]

    set_seq_count[s] += 1
    mult = 3
    if freq <= 2000:
      mult = -math.exp(-freq/300)*2.2+3

    data['width'] = mult*w_i

  add_edge_attributes(g, draw_mode)
  dot_file_args = {
      'mode' : 'w+',
      'prefix' : os.path.basename(output_name),
      'suffix' : '.dot',
      'delete' : not keep_tmp_files,
      'dir' : '.' if keep_tmp_files else None
  }
  with NamedTemporaryFile(**dot_file_args) as dot_fd:
    render_to_file(
        '/kstep.tmpl.dot',
        dot_fd,
        filters={ 'trimFn': lambda fn: os.path.splitext(os.path.basename(fn))[0] },
        context = {
          'color_scheme' : color_scheme,
          'legend_font_size' : determineLegendFontSize(g),
          'nodes' : g.nodes(data=True),
          'edges' : g.edges(data=True),
          'files' : [ (set_id,colors[set_id],set_seq_count[set_id]) for set_id in sorted(sample_names, key=lambda s: '' if s == 'shared' else s) ]
    })
    dot_fd.flush()

    try:
      subprocess.check_call(['neato', '-T', output_type, dot_fd.name, '-o', output_name + '.' + output_type])
    except OSError as ex:
      if ex.errno == 2:
        six.raise_from(
          OSError('Probabbly could not find neato. Check the graphviz installation.'),
          ex)
      else:
        six.reraise(*sys.exc_info())

def draw_kStep_main(args):
  g = nx.Graph()
  def sequences():
    freq_threshold = args.frequencycutoff

    seq_id=0
    for f_id,fn in enumerate(args.files, start=2):
      seq_count = 0
      set_id = first(os.path.basename(fn).split('.', 1))
      with open(fn) as in_fd:
        for seq in SeqIO.parse(in_fd, 'fasta'):
          try:
            freq = int(seq.description.rsplit('_',1)[1])
          except ValueError,IndexError:
            warnings.warn('Could not find frequency in sequence name')
            freq = 1

          if freq < freq_threshold:
            continue

          g.add_node(
            seq_id,
            sample=str(set_id),
            freq=freq)

          seq_count += 1
          seq_id += 1
          yield seq

  with msa.align(args.alignment_strategy, sequences()) as aln:
    draw_kstep(g,
      imap(lambda t: (t[2],(t[0],t[1])),
        chunked_hamming(aln, ignore_gaps=False, sort=True, key=itemgetter(2))),
      args.draw_mode,
      args.output,
      keep_tmp_files=args.keep_files,
      output_type=args.output_type)

def kstep_subcommand(subparsers):
  parser = subparsers.add_parser('kstep',
    description='Show k-step network for a given pair of files')
  parser.set_defaults(command=draw_kStep_main)
  parser.add_argument('files',
    nargs='+',
    help='List of files to be analyzed, order of files does not matter')
  parser.add_argument('-o', '--output',
    required=False, default='kstep',
    help='Your preferred output name (no file extension necessary)')
  parser.add_argument('-t', '--type',
    dest='output_type',
    default='png', choices=['dot', 'ps', 'svg', 'svgz', 'fig', 'mif', 'hpgl', 'pcl', 'pdf', 'png', 'gif', 'dia', 'imap', 'cmapx'],
    help='The output file type')
  parser.add_argument('-f', '--frequencycutoff',
    type=int, required=False, default=0,
    help='Minimum frequency sequence must have to appear on the figure, '
         'incompatible with -a')
  parser.add_argument('-a', '--dontalign',
    action='store_false', default=True,
    help='Pass this as an argument to skip file alignment; passing -a will '
    'override -f and will process all sequences on unaligned files')
  parser.add_argument('-k', '--keepfiles',
    action='store_true', default=False, dest='keep_files',
    help='Pass this as an argument to not remove the intermediate files') # and the dot file?
  parser.add_argument('-c', '--colorscheme',
    dest='graphviz_colorscheme', default='set28',
    help='The graphviz colorscheme used. There should be at least one more '
    'color than files. [default: %(default)s]')
  parser.add_argument('-d', '--drawmode',
    dest='draw_mode', default='dont', choices=['dont','red','white','allred','weights'],
    help='''Decide what to do with longer edges on the network. [default: %(default)s]
    #dont - remove edges from network, draw disconnected components if necessary - the disconnected components will not be placed in any specific way relative to one another by graphviz
    #allred - color longer edges red with their length set to 9 for the graphing software
    #red - only act on longer edges if they are between two different patients
    #white - all edges are drawn with their true weights, with longer links turned white
    #weights - draw all edges with labels (very busy chart), shorten longer edges''')
  parser.add_argument('-s', '--strategy',
      dest='alignment_strategy', default=msa.DEFAULT_STRATEGY_NAME,
      choices=msa.strategy_names(),
      help='The name of the strategy to use [default: %(default)s]')
