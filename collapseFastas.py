from Bio import SeqIO
from collections import defaultdict

def parse(files):
    seqs=defaultdict(int)
    for file in files:
        with open(file) as f:
            for record in SeqIO.parse(f,'fasta'):
                try:
                    freq=int(record.id.split('_')[-1])
                except:
                    newid=''
                    i=-1
                    char=record.id[i]
                    while char not in '0123456789':
                        i-=1
                        char=record.id[i]
                    while char in '0123456789':
                        i-=1
                        newid+=char
                        char=record.id[i]
                    try:
                        freq=int(newid[::-1])
                    except ValueError:
                        import sys
                        sys.exit(record.id)
                seqs[record.seq]+=freq
    return seqs
    
def main(files):
    seqs=parse(files)
    for seq, freq in sorted(seqs.iteritems(),reverse=True,key=lambda (k,v):(v,k)):
        print('>seq_'+str(freq))
        print(seq)

if __name__=='__main__':
    import argparse # possible arguments to add: delta, nIter
    parser = argparse.ArgumentParser(description='combine ur files, keep frequencies but otherwise butcher names')
    parser.add_argument('files',
        nargs='+',
        help="List  of  files  to  be  analyzed,  order  of  files  does  not  matter")
    args = parser.parse_args()
    main(args.files)
