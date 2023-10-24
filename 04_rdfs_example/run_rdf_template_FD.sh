#!/bin/bash

initFile="lk_8_init.xyz"
trajFile="lk_8_example_traj.dcd"

vmd "${initFile}" "${trajFile}" -dispdev text -e rdf_FD.tk
