#!/bin/sh

echo "$*"

mpiexec="$1"
mpi_np="$2"
refout="$3"
name="$4"
shift 4

tmpout="${name}.tmp"
echo PWD: $PWD
echo TMPOUT: $tmpout
echo REFOUT: $refout
"${mpiexec}" -n ${mpi_np} "./${name}" "$@" > "${tmpout}" || exit 1
diff -u ${refout} ${tmpout} || exit 1
rm "${tmpout}"
