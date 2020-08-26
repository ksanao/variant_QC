#!/usr/bin/env bash
# Author: Oksana Riba Grognuz
# Creation Date: 2020.08.21

ERROR_ARGUMENTS=65

USAGE="\nUsage: \n\
`basename $0` [-n container_name] run\n\
 \tdefault container_name=variant_qc\n"

# Print usage: scriptname -options, if script invoked with no command-line args
if [ $# -eq 0 ]  
then
	echo -e $USAGE
	exit $ERROR_ARGUMENTS
fi  

# Options followed by : expect an argument 
while getopts "n:" Option
do
  case $Option in
    n ) container_name=$OPTARG;;
    h ) echo -e $USAGE; exit 1;;
    #if no argument is supplied to the option requiring it, it will fall into default
    * ) echo "Unimplemented option chosen.";  exit $ERROR_ARGUMENTS;;   # DEFAULT
  esac
done

# Decrements the argument pointer so it points to next argument.
shift $(($OPTIND - 1))

#### DEFAULT VALUES
container_name_DEFAULT="variant_qc"

#### SETTING UNDEFINED OPTIONS TO DEFAULTS
: ${container_name=$container_name_DEFAULT}

#chmod +x variant_QC/docker-entrypoint.sh
docker run -p 8888:8888 --name $container_name --env=/root/.profile -ti -v $PWD/variant_QC/:/home/ variant_qc

exit 0

