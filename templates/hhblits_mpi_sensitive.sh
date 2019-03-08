# Params
# NTASKS
# CPU_PER_TASK
# INPUT
# OUTPUT
# DB

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
      -n 3 \
      -loc \
      -realign_old_hits \
      -M a2m \
      -mact 0.4 \
      -v 1 \
      -d "${DB}/db"
