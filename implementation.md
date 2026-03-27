# Phase 87 — NetworkX Bipartite Graph (v1.0)

## Goal
Build a bipartite graph connecting compounds (from compounds.csv) to the KRAS target, visualize with matplotlib.

## CLI
```bash
PYTHONUTF8=1 python main.py --compounds data/compounds.csv
```

## Outputs
- `output/bipartite_graph.png` — bipartite graph visualization
- Console summary of graph statistics

## Logic
1. Load compounds.csv (compound_name, smiles, pic50)
2. Create bipartite graph: compound nodes (left) + KRAS target node (right)
3. Edge weight = pIC50 value
4. Color compounds by scaffold family (benz/naph/ind/quin/pyr/bzim/other)
5. Layout and save figure

## Key Concepts
- NetworkX bipartite graph construction
- Scaffold family color palette
- Graph visualization with matplotlib

## Verification Checklist
- [ ] `--help` works
- [ ] Compounds loaded correctly
- [ ] Graph has correct node/edge counts
- [ ] PNG saved with bipartite layout
- [ ] Scaffold family colors applied

## Risks
- Large number of compounds may make graph cluttered (45 is manageable)
