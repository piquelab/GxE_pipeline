#!/bin/bash
set -e
set -v

mkdir -p ../../expression/fpkms/
mkdir -p ../../expression/DESeq2/

cp util_expression_runAll.sh ../../expression/fpkms/
cp util_combine.sh ../../expression/fpkms/
cp expr_* ../../expression/fpkms/

cp util_DEG_analysis.sh ../../expression/DESeq2/
cp DESeq2_* ../../expression/DESeq2/

cp treatmentKey.txt ../../expression/
cp plateKey.txt ../../expression/
