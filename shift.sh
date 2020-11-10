
#!/bin/bash
# shft.sh: Using 'shift' to step through all the positional parameters.

#  Name this script something like shft.sh,
#+ and invoke it with some parameters.
#+ For example:
#             sh shft.sh a b c def 83 barndoor

until [ -z "$1" ]  # Until all parameters used up . . .
do
  echo -n "$1 "
  shift
  echo
done


