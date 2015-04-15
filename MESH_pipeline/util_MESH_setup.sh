#!/bin/bash
set -v
set -e

## copy the script for assembling the master table
mkdir -p ../../jointGenotyping/MESH_results_masterTable/data_logFC/
cp MESH_QuASAR_master* ../../jointGenotyping/MESH_results_masterTable/data_logFC

