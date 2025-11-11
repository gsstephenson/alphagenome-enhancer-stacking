#!/bin/bash
# Logic Gates Experiment - Master Pipeline
# Run all steps sequentially

set -e  # Exit on error

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
cd "$SCRIPT_DIR"

echo "================================"
echo "LOGIC GATES EXPERIMENT PIPELINE"
echo "================================"
echo ""

# Step 1: Build constructs
echo "[1/4] Building logic gate constructs..."
python build_logic_gate_constructs.py
echo "✓ Step 1 complete"
echo ""

# Step 2: Run predictions
echo "[2/4] Running AlphaGenome predictions..."
echo "NOTE: This may take several hours depending on your AlphaGenome setup"
read -p "Press Enter to continue or Ctrl+C to abort..."
python run_logic_gate_predictions.py
echo "✓ Step 2 complete"
echo ""

# Step 3: Analyze results
echo "[3/4] Analyzing logic gate patterns..."
python analyze_logic_gates.py
echo "✓ Step 3 complete"
echo ""

# Step 4: Generate figures
echo "[4/4] Creating visualization figures..."
python create_logic_gate_figures.py
echo "✓ Step 4 complete"
echo ""

echo "================================"
echo "✓ PIPELINE COMPLETE!"
echo "================================"
echo ""
echo "Results saved to:"
echo "  - Constructs: sequences/"
echo "  - Predictions: alphagenome_outputs/"
echo "  - Analysis: logic_gate_analysis.json"
echo "  - Figures: results/figures/"
echo ""
echo "Next steps:"
echo "  1. Review results/figures/summary_dashboard.png"
echo "  2. Check logic_gate_analysis.json for detailed metrics"
echo "  3. Write manuscript based on findings!"
