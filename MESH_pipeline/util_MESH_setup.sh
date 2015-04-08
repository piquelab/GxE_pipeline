#!/bin/bash
set -v
set -e

## create analysis directory & copy scripts into the directory for MESH analysis
mkdir -p ../../jointGenotyping/MESH_results_masterTable/analysis/
cp ./MESH/* ../../jointGenotyping/MESH_results_masterTable/analysis/

## copy the script for assembling the master table
mkdir -p ../../jointGenotyping/MESH_results_masterTable/data_logFC/
cp MESH_QuASAR_master* ../../jointGenotyping/MESH_results_masterTable/data_logFC

