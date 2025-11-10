# GitHub Copilot Instructions for mocbapy

## Project Overview

**mocbapy** (Multi-Objective Constraint-Based Analysis in Python) is a computational biology package that implements Multi-Objective Flux Balance Analysis (MO-FBA) and Multi-Objective Flux Variability Analysis (MO-FVA) for metabolic modeling. The implementation is based on the methodology described in [Budinich (2017)](https://doi.org/10.1371/journal.pone.0171744).

### Key Purpose
mocbapy enables the analysis of microbial ecosystems by creating multi-compartment models where multiple organisms (represented as COBRA models) interact through a shared metabolite pool. It computes Pareto-optimal solutions for competing biological objectives using vector linear programming.

## Architecture and Dependencies

### Core Dependencies
- **benpy**: Vector Linear Programming (VLP) solver that implements the bensolve algorithm. This is the mathematical engine for multi-objective optimization.
- **COBRApy (cobra)**: Constraint-Based Reconstruction and Analysis toolbox for Python. Provides metabolic model representation and basic FBA capabilities.
- **optlang**: Abstraction layer for mathematical optimization solvers (GLPK, Gurobi, CPLEX, etc.)
- **scipy**: Sparse matrix operations
- **numpy**: Numerical computations
- **pandas**: Data structure for FVA results
- **tqdm**: Progress bars for long-running operations

### Relationship with benpy
- mocbapy **depends on** benpy; it is NOT a standalone solver
- benpy provides the `vlpProblem` class and `solve()` function for multi-objective optimization
- mocbapy transforms ecosystem models into VLP problems that benpy can solve
- The conversion happens in `EcosystemModel.to_vlp()` method

### COBRApy Integration
- Individual organisms are represented as COBRApy `Model` objects
- mocbapy uses COBRA's:
  - Model structure (reactions, metabolites, bounds)
  - Stoichiometric matrix creation (`create_stoichiometric_matrix`)
  - Exchange reactions identification (`model.exchanges`)
  - Solver interface abstractions (via optlang)
- Does NOT use COBRA's optimization functions directly (uses benpy instead)

## Project Structure

```
src/mocbapy/
├── __init__.py           # Package initialization
├── EcosystemModel.py     # Core class for multi-compartment models
├── analysis.py           # MO-FBA and MO-FVA implementations
├── utilities.py          # Helper functions (solver selection, model building)
├── draw.py              # Visualization wrappers for benpy.draw
├── test2D.py            # 2D example/test script
└── test3D.py            # 3D example/test script

tests/
└── TestFBA.py           # FBA testing with multiple SBML models
```

## Core Concepts

### EcosystemModel
The central class that represents a multi-organism metabolic ecosystem:
- **model_array**: List of COBRApy models representing different organisms
- **metabolic_dict**: Dictionary mapping `(metabolite_id, model)` tuples to shared metabolite names in the pool
- **medium**: Optional dictionary defining nutrient availability with `{metabolite: (lower_bound, upper_bound)}`
- **pool**: Shared metabolite compartment where organisms exchange metabolites

### MO-FBA (Multi-Objective Flux Balance Analysis)
- Simultaneously optimizes growth of multiple organisms
- Returns a Pareto front of non-dominated solutions
- Each organism's biomass objective is treated as a separate objective
- Uses benpy's VLP solver to compute the entire Pareto front

### MO-FVA (Multi-Objective Flux Variability Analysis)
- Determines flux ranges near Pareto-optimal solutions
- Given a point on/near the Pareto front, constrains objectives and explores flux space
- Uses single-objective optimization (via optlang) iteratively
- `alpha` parameter controls proximity to Pareto front (0 < alpha ≤ 1)

## Coding Standards

### Language and Style
- **Pure Python 3.9+**: No Cython, no compiled extensions
- Follow PEP 8 style guidelines
- Use descriptive variable names (e.g., `ecosystem_model`, not `em`)
- UTF-8 encoding with `# -*- coding: utf-8 -*-` headers in older files

### Code Patterns to Follow

#### Import Organization
```python
# Standard library first
from collections import OrderedDict, defaultdict
from warnings import warn

# Third-party packages
from benpy import vlpProblem
import cobra
import numpy as np
from scipy.sparse import lil_matrix

# Local imports
from .utilities import choose_optlang_interfase, build_base_opt_model
```

#### Sparse Matrix Usage
Always use sparse matrices (scipy.sparse) for stoichiometric matrices:
```python
from scipy.sparse import lil_matrix, block_diag, eye

# Use lil_matrix for construction
array_form = block_diag([create_stoichiometric_matrix(model, array_type="lil")
                         for model in self.models], format="lil")

# Access with standard indexing
self.Ssigma[met_idx, rxn_idx] = -coeff
```

#### Naming Conventions
- **Variables**: `snake_case` (e.g., `ecosystem_model`, `pool_ex_rxns`)
- **Functions/Methods**: `snake_case` (e.g., `mo_fba`, `to_vlp`)
- **Classes**: `PascalCase` (e.g., `EcosystemModel`)
- **Private methods**: Leading underscore (e.g., `_construct_ecosystem_pool`)
- **Organism identifiers**: Appended with `:{model.id}` (e.g., `"BIOMASS:ecoli_1"`)
- **Pool metabolites**: Appended with `:pool` (e.g., `"glc__D_e:pool"`)

#### Error Handling
Use warnings for non-critical issues, exceptions for errors:
```python
from warnings import warn

# Non-critical issues
if solver is None:
    warn("Warning: No solver selected. Available solvers: {}".format(str(avail)))

# Critical errors
if name in self.pool_ex_rxns:
    raise RuntimeError("Reaction already included")
```

#### Progress Indicators
Use tqdm for long-running operations:
```python
from tqdm import tqdm

for rxn in tqdm(reactions):
    # Long operation
```

### Type Expectations
- Model arrays: List of COBRApy `Model` objects
- Metabolic dictionaries: `dict` with `(str, Model)` tuples as keys
- Bounds: Two-element tuples `(lower_bound, upper_bound)`
- Stoichiometric matrices: scipy sparse matrices (lil_matrix)
- Objective dictionaries: `OrderedDict` with Model objects as keys

## FBA and MO-FBA Guidelines

### Standard FBA (Single Organism)
When working with individual COBRA models:
```python
import cobra

model = cobra.io.read_sbml_model('path/to/model.xml')
solution = model.optimize()  # Uses COBRApy's optimizer
```

### MO-FBA (Multiple Organisms)
When working with ecosystem models:
```python
from mocbapy import EcosystemModel, mo_fba, utilities

# Create ecosystem
models = [model1, model2, model3]
common_mets = utilities.get_common_mets(models)
ecosystem = EcosystemModel(models, common_mets, medium)

# Solve with benpy
solution = mo_fba(ecosystem, options=bensolve_options)

# Access Pareto front
primal_vertices = solution.Primal.vertex_value
```

### Critical FBA Concepts
1. **Stoichiometric Matrix (S)**: Rows = metabolites, columns = reactions, values = stoichiometric coefficients
2. **Mass Balance**: `Sv = 0` (steady-state assumption)
3. **Flux Bounds**: `lb ≤ v ≤ ub` for each reaction
4. **Objective Function**: Typically biomass reaction (growth rate)
5. **Exchange Reactions**: Interface between model and environment (e.g., `EX_glc__D_e`)

### Optimization Solver Hierarchy
1. Commercial solvers (if available): Gurobi > CPLEX
2. Open-source: GLPK (via optlang)
3. Fallback: scipy (slow, not recommended for large models)

## Testing Standards

### Test Organization
- Tests go in `tests/` directory
- Use descriptive test file names (e.g., `TestFBA.py`)
- Include example/test scripts in `src/mocbapy/` (e.g., `test2D.py`, `test3D.py`)

### Test Models
Use COBRA's built-in test models when possible:
```python
import cobra.test
model = cobra.test.create_test_model("ecoli")  # or "textbook"
```

For ecosystem tests, create multiple instances:
```python
models = []
for i in range(n_organisms):
    model = cobra.test.create_test_model("ecoli")
    model.id = f"{model.id}_{i+1}"  # Unique IDs required
    models.append(model)
```

### Test Coverage Areas
1. **Model Construction**: EcosystemModel initialization
2. **Pool Creation**: Shared metabolite space
3. **Constraint Handling**: Bounds, medium constraints
4. **MO-FBA Execution**: Pareto front computation
5. **MO-FVA Execution**: Flux ranges near Pareto solutions
6. **Edge Cases**: Empty models, missing dictionaries, single organism

### Validation Patterns
```python
# Check model dimensions
assert len(ecosystem.sysreactions) == sum(len(m.reactions) for m in models) + len(pool_reactions)
assert len(ecosystem.sysmetabolites) == sum(len(m.metabolites) for m in models) + len(pool_metabolites)

# Check stoichiometric matrix shape
m, n = ecosystem.Ssigma.shape
assert m == len(ecosystem.sysmetabolites)
assert n == len(ecosystem.sysreactions)

# Verify solution structure
assert hasattr(solution, 'Primal')
assert hasattr(solution, 'Dual')
```

### Performance Considerations
- FBA problems scale with number of reactions (typically O(n³))
- MO-FBA problems scale with number of objectives and reactions
- Use solver='gurobi' or 'cplex' for large models (>1000 reactions)
- Cache optimization models when performing repeated FVA

## Common Pitfalls and Best Practices

### Pitfalls to Avoid
1. **Non-unique model IDs**: Always set unique `model.id` for each organism
2. **Inconsistent metabolite naming**: Use metabolic_dict to map metabolites correctly
3. **Missing solvers**: Check solver availability with `optlang.available_solvers`
4. **Dense matrix operations**: Always use sparse matrices for stoichiometry
5. **Unbounded pools**: Set reasonable bounds in medium dictionary
6. **Ignoring alpha in MO-FVA**: alpha=1 means exactly on Pareto front (may be infeasible)

### Best Practices
1. **Validate input models**: Check for orphan reactions/metabolites before building ecosystem
2. **Use progress indicators**: Wrap long loops with `tqdm()`
3. **Set bensolve options**: Control verbosity with `options['message_level']`
4. **Normalize Pareto plots**: Use max values for better visualization
5. **Document metabolic dictionaries**: Comment mapping logic for shared metabolites

## Common Tasks and Patterns

### Creating an Ecosystem Model
```python
from mocbapy import EcosystemModel, utilities

# Load individual models
models = [cobra.io.read_sbml_model(path) for path in model_paths]

# Automatic common metabolite detection
common_mets = utilities.get_common_mets(models)

# Define medium (optional)
medium = {
    'glc__D_e': (-10, 1000),  # Glucose limited
    'o2_e': (-20, 1000),      # Oxygen
}

# Create ecosystem
ecosystem = EcosystemModel(models, common_mets, medium)
```

### Running MO-FBA and Extracting Results
```python
from mocbapy import mo_fba, bensolve_default_options
import numpy as np

# Configure bensolve
options = bensolve_default_options()
options['message_level'] = 0  # Suppress output

# Solve
solution = mo_fba(ecosystem, options=options)

# Extract Pareto front vertices (extreme points)
ext_points = np.asarray([
    v for i, v in enumerate(solution.Primal.vertex_value)
    if solution.Primal.vertex_type[i] == 1
])

# Get maximum growth for each organism
max_growth = np.amax(ext_points, axis=0)
```

### Running MO-FVA
```python
from mocbapy import mo_fva

# Define point on/near Pareto front
fba_result = {
    'BIOMASS:organism1': 0.5,
    'BIOMASS:organism2': 0.3,
}

# Run FVA on exchange reactions
exchange_rxns = [r for r in ecosystem.sysreactions if ':pool' in r]
fva_result = mo_fva(ecosystem, fba=fba_result, reactions=exchange_rxns, alpha=0.9)

# fva_result is a pandas DataFrame with 'minimum' and 'maximum' columns
```

### Adding Manual Pool Reactions
```python
# Add a cross-feeding reaction
reaction_dict = {
    'acetate:pool': -1,  # Consumed from pool
    'co2:pool': 2,       # Produced to pool
}
ecosystem.add_pool_reaction(
    name='Custom_CrossFeed:pool',
    reaction=reaction_dict,
    lb=0,    # Irreversible
    ub=100
)
```

## Development Workflow

1. **Load Models**: Use COBRA to load/create individual organism models
2. **Define Connections**: Create metabolic_dict for shared metabolites
3. **Build Ecosystem**: Instantiate EcosystemModel
4. **Set Constraints**: Adjust bounds (medium, exchange reactions)
5. **Solve MO-FBA**: Get Pareto front with benpy
6. **Analyze Results**: Extract extreme points, visualize
7. **Optional FVA**: Explore flux ranges near solutions

## Documentation Expectations

- Docstrings for all public functions and classes
- Explain multi-objective optimization concepts when relevant
- Reference Budinich (2017) paper for methodology details
- Include example usage in docstrings for complex functions
- Document parameter units and expected types

## Version Compatibility

- Python: 3.9 - 3.12
- COBRApy: Compatible with current releases (0.25+)
- benpy: As specified in requirements
- optlang: As dependency of COBRApy

## Additional Resources

- **Original Paper**: Budinich M, et al. (2017) "A multi-objective constraint-based approach for modeling genome-scale microbial ecosystems" PLOS ONE
- **COBRApy Documentation**: https://cobrapy.readthedocs.io/
- **benpy**: Multi-objective optimization via bensolve
- **FBA Review**: Orth JD, et al. (2010) "What is flux balance analysis?" Nature Biotechnology
