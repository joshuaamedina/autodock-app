# Allow over-ride
if [ -z "${CONTAINER_IMAGE}" ]
then
    version=$(cat ./_util/VERSION)
    CONTAINER_IMAGE="index.docker.io/library/ubuntu:bionic"
fi
. lib/container_exec.sh

# Write an excution command below that will run a script or binary inside the 
# requested container, assuming that the current working directory is 
# mounted in the container as its WORKDIR. In place of 'docker run' 
# use 'container_exec' which will handle setup of the container on 
# a variety of host environments. 
#
# Here is a template:
#
# container_exec ${CONTAINER_IMAGE} COMMAND OPTS INPUTS
#
# Here is an example of counting words in local file 'poems.txt',
# outputting to a file 'wc_out.txt'
#
# container_exec ${CONTAINER_IMAGE} wc poems.txt > wc_out.txt
#

# set -x

# set +x


singularity pull vina_1.2.3.sif docker://joshuaamedina2000/autodock-vina:1.2.3
flex=${flexible_sidechains}
if
     [ -z "${flexible_sidechains}" ]
then
     flex='Empty'
fi
    
MV2_ENABLE_AFFINITY=0 ibrun -np 32 singularity exec vina_1.2.3.sif python3 autodock.py -r ${receptor} -c "${center_x},${center_y},${center_z}" -s "${size_x},${size_y},${size_z}" -m ${forcefield} -d ${docking} -ll ${library} -n ${top_n_scores} -f $flex