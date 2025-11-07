#!/bin/bash

# Root directory
ROOT="/mnt/work_1/gest9386/CU_Boulder/rotations/LAYER/AlphaGenome_EnhancerStacking"

echo "📁 Creating project directories..."
mkdir -p $ROOT/{data/{raw,processed},filler,logs}
mkdir -p $ROOT/sequences/{enhancers,promoters,constructs}
mkdir -p $ROOT/alphagenome/{inputs,outputs}
mkdir -p $ROOT/analysis/{code,results}

# Step 1: Download HS2 enhancer sequence (chr11:5290000–5291000, GRCh38)
echo "⬇️  Downloading HS2 enhancer sequence..."
curl -s "http://genome.ucsc.edu/cgi-bin/das/hg38/dna?segment=chr11:5290000,5291000" \
  -H "Accept: text/plain" \
  -o $ROOT/sequences/enhancers/HS2_enhancer.fa

# Step 2: Download HBG1 promoter sequence (chr11:5273600–5273900, GRCh38)
echo "⬇️  Downloading HBG1 promoter sequence..."
curl -s "http://genome.ucsc.edu/cgi-bin/das/hg38/dna?segment=chr11:5273600,5273900" \
  -H "Accept: text/plain" \
  -o $ROOT/sequences/promoters/HBG1_promoter.fa

# Step 3: Generate 1M bp A/T-rich filler sequence
echo "🧬 Generating 1M A/T-rich filler sequence..."
python3 - <<EOF
import random
bases = ["A", "T", "G", "C"]
weights = [0.4, 0.4, 0.1, 0.1]  # A/T-rich profile
seq = "".join(random.choices(bases, weights=weights, k=1_000_000))
with open("$ROOT/filler/1M_filler.txt", "w") as f:
    f.write(seq)
EOF

# Step 4: Confirm final structure
echo "📂 Final directory tree:"
tree -L 3 $ROOT

echo "✅ Setup complete!"
