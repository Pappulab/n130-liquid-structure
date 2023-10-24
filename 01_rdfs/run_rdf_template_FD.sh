#!/bin/bash

initFile="##TOPODIR##/##SYSNAME##_init.xyz"
trajFile="##RUNSDIR##/##SYSNAME##/Run_##REPSNUM##/##SYSNAME##_##REPSNUM##_nvt.dcd"

vmd "${initFile}" "${trajFile}" -dispdev text -e rdf_FD.tk
