@echo off
docker --context desktop-linux run --rm -v "%CD%:/work" -w /work quay.io/biocontainers/miniprot:0.18--h577a1d6_0 miniprot %*
