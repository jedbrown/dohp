#!/bin/sh

echo "$*"

mpiexec="$1"
mpi_np="$2"
refout="$3"
name="$4"
shift 4

tmpout="${name}.tmp"
echo :: "${mpiexec}" -n ${mpi_np} "./${name}" "$@" || exit 1
"${mpiexec}" -n ${mpi_np} "./${name}" "$@" > "${tmpout}" || exit 1
diff -u ${refout} ${tmpout} || exit 1
rm "${tmpout}"
