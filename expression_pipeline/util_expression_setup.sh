#!/bin/bash
set -e
set -v

mkdir -p ../../expression/fpkms/
mkdir -p ../../expression/DESeq2/

cp expr_* ../../expression/fpkms/
cp DESeq2_* ../../expression/DESeq2/

cp treatmentKey.txt ../../expression/
cp plateKey.txt ../../expression/
