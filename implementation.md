# Phase 87 — NetworkX Bipartite Graph (v1.1)

## Goal
Build a bipartite graph connecting compounds (from compounds.csv) to the KRAS target, visualize with matplotlib.

## CLI
```bash
PYTHONUTF8=1 python main.py --compounds data/compounds.csv
```

## Outputs
- `output/bipartite_graph.png` — bipartite graph with scaffold family coloring
- Console summary of graph stats and family breakdown

## Logic
1. Load compounds.csv (45 compounds: compound_name, smiles, pic50)
2. Create bipartite graph: compound nodes (bipartite=0) + KRAS (bipartite=1)
3. Edge weight = pIC50; edge alpha proportional to pIC50
4. Color compounds by scaffold family using FAMILY_COLORS palette
5. Bipartite layout: compounds left, KRAS right

## Key Concepts
- NetworkX bipartite graph construction and validation
- Scaffold family color palette (established convention)
- Graph visualization with custom layout

## Verification Checklist
- [x] `--help` works
- [x] 45 compounds loaded correctly
- [x] Graph: 45 compounds, 1 target, 45 edges, is_bipartite=True
- [x] PNG saved with bipartite layout
- [x] Scaffold family colors applied (benz:12, naph:7, ind:7, quin:7, pyr:6, bzim:6)

## Results
- 45 compounds connected to KRAS via pIC50-weighted edges
- 6 scaffold families represented, benz most common (12)
- Graph confirmed bipartite via networkx.algorithms.bipartite.is_bipartite()

## Deviations
- None

## Risks
- Large compound sets would need different layout strategy
