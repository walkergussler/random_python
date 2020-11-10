import argparse

servers = [ "ApaServer", "BananServer", "GulServer", "SolServer", "RymdServer",
            "SkeppServer", "HavsServer", "PiratServer", "SvartServer", "NattServer", "SovServer" ]

parser = argparse.ArgumentParser(description="A program to update components on servers.")
group = parser.add_mutually_exclusive_group()
group.add_argument('-l', '--list', dest="update", action='store_false', default=False, help='list server components')
group.add_argument('-u', '--updatepom', dest="update", action='store_true', help='update server components')
parser.add_argument('-o', '--only', choices=servers, help='Space separated list of case sensitive server names to process.  Allowed values are '+', '.join(servers), metavar='')
parser.add_argument('-s', '--skip', choices=servers, help='Space separated list of case sensitive server names to exclude from processing.  Allowed values are '+', '.join(servers), metavar='')
args = parser.parse_args()

