#!/bin/bash

initFile="rpL5_init.xyz"
trajFile="rpL5_example_traj.dcd"

vmd "${initFile}" "${trajFile}" -dispdev text -e rdf_PD.tk
