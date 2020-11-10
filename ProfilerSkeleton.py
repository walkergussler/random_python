import pprofile

prof=pprofile.Profile()
with prof():

#whole program

prof.print_stats()


