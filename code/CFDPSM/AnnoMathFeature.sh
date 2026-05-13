#!/bin/bash

set -e

base_path="/MathFeature/" # where you install the mathfeature
outdirPath="$base_path/"
ref_inputfile="$outdirPath/ref_sequences.fasta"
alt_inputfile="$outdirPath/alt_sequences.fasta"

echo "[INFO] Processing SE..."
outfilePath3_1="$outdirPath/SE_alt.csv"
python3.7 $base_path/methods/EntropyClass.py -i "$alt_inputfile" -o "$outfilePath3_1" -l DNA -k 5 -e Shannon
outfilePath3_2="$outdirPath/SE_ref.csv"
python3.7 $base_path/methods/EntropyClass.py -i "$ref_inputfile" -o "$outfilePath3_2" -l DNA -k 5 -e Shannon


echo "[INFO] Processing Chaos Game Representation(CCGR)..."
outfilePath2_1="$outdirPath/CCGR.csv"
python3.7 $base_path/methods/ChaosGameTheory_sigma.py -o "$outfilePath2_1" -r 1 -d "${alt_inputfile}:0"

python3.7 Filter_feature_MathFeature.py --i1 "$outfilePath3_2" --i2 "$outfilePath3_1" --i3 "$outfilePath2_1" --output "$outdirPath/mathfeature.csv"