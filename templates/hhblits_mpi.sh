# Params
# NTASKS
# CPU_PER_TASK
# INPUT
# OUTPUT
# DB
# NITER

mpirun -np ${NTASKS} -cpus-per-proc ${CPU_PER_TASK} \
  ffindex_apply_mpi \
    ${INPUT}.ff{data,index} \
    -i "${OUTPUT}.ffindex" \
    -d "${OUTPUT}.ffdata" \
    -- \
    hhblits \
      -i stdin \
      -o stdout \
      -cpu ${CPU_PER_TASK} \
      -n ${NITER} \
      -loc \
      -M a2m \
      -mact 0.4 \
      -v 1 \
      -d "${DB}/db"
