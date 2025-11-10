# Agent Coordination Guide

This document helps you choose the right specialized agent for different tasks in the mocbapy project.

## Available Agents

### 1. Python Developer (`dev_base`)
**Use for**: General Python development, new features, API design
- Implementing new Python features
- Refactoring pure Python code
- COBRApy and benpy integration
- Creating or modifying example scripts
- General bug fixes in Python code
- Extending EcosystemModel or analysis functions

### 2. Bioinformatics & Metabolic Modeling Expert (`bioinformatics-agent`)
**Use for**: Biological validation and interpretation
- Validating COBRA models for biological correctness
- Designing multi-organism ecosystems
- Interpreting MO-FBA/MO-FVA results biologically
- Setting up realistic medium constraints
- Identifying cross-feeding patterns
- Validating growth rates and flux ranges
- Understanding Pareto fronts in biological context
- Connecting results to biological literature

### 3. Documentation Expert (`docagent`)
**Use for**: All documentation-related tasks
- Writing or updating README
- Creating usage examples with biological context
- Maintaining CHANGELOG
- Writing API docstrings
- Creating installation guides
- Explaining biological concepts
- Documenting metabolic modeling workflows

### 4. Testing Expert (`testing-agent`)
**Use for**: Testing and quality assurance
- Creating test suites (unit, integration)
- Writing tests for new features
- Validating numerical accuracy
- Testing biological model validity
- Creating COBRA model test fixtures
- Testing edge cases and error handling
- Setting up test infrastructure
- Improving test coverage

### 5. Security & Code Quality Expert (`security-agent`)
**Use for**: Security audits and code quality
- Finding and fixing input validation issues
- Code security reviews
- Dependency vulnerability scanning
- Static analysis and linting
- Improving error handling
- Code quality improvements
- Type safety and type hints

## Task-to-Agent Mapping

### Common Scenarios

| Task | Primary Agent | Supporting Agents |
|------|---------------|-------------------|
| Add new MO-FBA feature | dev_base | testing-agent, docagent |
| Validate biological model | bioinformatics-agent | - |
| Interpret Pareto front | bioinformatics-agent | docagent |
| Update README | docagent | - |
| Add unit tests | testing-agent | - |
| Design ecosystem | bioinformatics-agent | dev_base |
| Fix bug in Python code | dev_base | testing-agent |
| Security review | security-agent | - |
| Add API documentation | docagent | dev_base |
| Validate model inputs | security-agent | bioinformatics-agent |
| Create example script | dev_base | bioinformatics-agent, docagent |
| Improve code quality | security-agent | dev_base |

### Multi-Agent Workflows

#### New Feature Development
1. **dev_base**: Implement the feature
2. **bioinformatics-agent**: Validate biological correctness
3. **testing-agent**: Create comprehensive tests
4. **security-agent**: Review for security issues
5. **docagent**: Document the new feature

#### Biological Model Analysis
1. **bioinformatics-agent**: Design and validate ecosystem model
2. **dev_base**: Implement any custom model building code
3. **testing-agent**: Create tests for model construction
4. **docagent**: Document the biological rationale

#### Bug Fix Workflow
1. **testing-agent**: Create failing test that reproduces bug
2. **dev_base**: Fix the bug
3. **bioinformatics-agent**: Validate biological correctness (if applicable)
4. **security-agent**: Check if fix introduces vulnerabilities

#### Documentation Update
1. **docagent**: Write/update documentation
2. **bioinformatics-agent**: Validate biological accuracy
3. **dev_base**: Ensure code examples work
4. **testing-agent**: Convert examples to tests if appropriate

#### Ecosystem Modeling Task
1. **bioinformatics-agent**: Design ecosystem structure and constraints
2. **dev_base**: Implement ecosystem construction code
3. **testing-agent**: Validate with test COBRA models
4. **docagent**: Document the ecosystem design and interpretation

## Detailed Task Breakdowns

### Adding a New Analysis Function

**Example**: Add a new function to compute metabolite exchange patterns

1. **dev_base** (Implementation)
   - Write the core Python function
   - Integrate with EcosystemModel class
   - Handle sparse matrices efficiently
   - Add error handling

2. **bioinformatics-agent** (Validation)
   - Ensure biological correctness of algorithm
   - Validate with known metabolic patterns
   - Check against biological literature
   - Verify results make biological sense

3. **testing-agent** (Testing)
   - Create unit tests with COBRA test models
   - Test edge cases (empty models, single organism)
   - Validate numerical accuracy
   - Test integration with mo_fba/mo_fva

4. **security-agent** (Security Review)
   - Validate input parameters
   - Check for numerical stability
   - Ensure proper error handling
   - Review for potential vulnerabilities

5. **docagent** (Documentation)
   - Write comprehensive docstring
   - Add usage examples
   - Explain biological interpretation
   - Update README if necessary

### Validating an Ecosystem Model

**Example**: User creates an ecosystem but gets unexpected results

1. **bioinformatics-agent** (Primary)
   - Check model structure (orphan metabolites/reactions)
   - Validate stoichiometry and mass balance
   - Review medium constraints for realism
   - Interpret results in biological context
   - Suggest improvements

2. **dev_base** (Supporting)
   - Check if code is using mocbapy correctly
   - Verify COBRApy integration
   - Ensure proper use of exchange reactions

3. **testing-agent** (Validation)
   - Create test case to reproduce issue
   - Verify with simpler test models
   - Validate against expected behavior

### Creating Example Scripts

**Example**: Add a new example for three-species ecosystem

1. **bioinformatics-agent** (Design)
   - Choose biologically relevant organisms
   - Design realistic medium
   - Identify shared metabolites
   - Define biological question being addressed

2. **dev_base** (Implementation)
   - Write the example script
   - Use cobra.test models or load SBML
   - Implement model construction
   - Add plotting/visualization

3. **docagent** (Documentation)
   - Add detailed comments explaining biology
   - Explain each step clearly
   - Interpret results
   - Add to README or examples section

4. **testing-agent** (Validation)
   - Ensure example runs without errors
   - Verify it produces expected output
   - Consider converting to integration test

## Agent Coordination Best Practices

### 1. Start with the Right Agent
- For biological questions → **bioinformatics-agent**
- For code implementation → **dev_base**
- For testing → **testing-agent**
- For documentation → **docagent**
- For security/quality → **security-agent**

### 2. Chain Agents Sequentially
- Complete one agent's work before moving to next
- Pass context between agents clearly
- Validate each agent's output before proceeding

### 3. Use Multiple Agents for Complex Tasks
- Don't try to do everything with one agent
- Leverage specialized expertise
- Coordinate handoffs between agents

### 4. Prioritize by Domain
- **Biology-first**: For metabolic modeling tasks, consult bioinformatics-agent early
- **Code-first**: For technical implementation, start with dev_base
- **Testing-throughout**: Include testing-agent early in development

### 5. Document Agent Decisions
- Explain why you chose certain agents
- Note any conflicts between agent recommendations
- Document resolution of conflicts

## Common Agent Combinations

### Biological Analysis Workflow
```
bioinformatics-agent (design) 
  → dev_base (implement) 
  → testing-agent (validate) 
  → docagent (document)
```

### Technical Bug Fix
```
testing-agent (reproduce bug) 
  → dev_base (fix) 
  → security-agent (review) 
  → bioinformatics-agent (if biological validation needed)
```

### Documentation Enhancement
```
docagent (write) 
  → bioinformatics-agent (validate biology) 
  → dev_base (verify code examples)
```

### New Feature Development
```
bioinformatics-agent (biological requirements)
  → dev_base (implement)
  → security-agent (security review)
  → testing-agent (test)
  → docagent (document)
```

## When to Use Which Agent

### Use **bioinformatics-agent** when:
- Designing ecosystems
- Validating biological correctness
- Interpreting results
- Setting medium constraints
- Understanding cross-feeding
- Questioning growth rates
- Explaining Pareto fronts
- Connecting to biological literature

### Use **dev_base** when:
- Writing Python code
- Integrating COBRApy/benpy
- Refactoring
- Fixing bugs
- Adding features
- Working with sparse matrices
- Building models

### Use **testing-agent** when:
- Creating tests
- Debugging failures
- Validating numerical accuracy
- Testing edge cases
- Improving coverage
- Creating fixtures

### Use **docagent** when:
- Writing documentation
- Creating examples
- Explaining concepts
- Maintaining README
- Writing docstrings
- Creating tutorials

### Use **security-agent** when:
- Reviewing code quality
- Validating inputs
- Checking dependencies
- Improving error handling
- Running linters
- Security audits

## Special Considerations for mocbapy

### Biology Matters
- **Always** consider biological validity
- Consult **bioinformatics-agent** for model design decisions
- Validate that mathematical solutions make biological sense

### Pure Python Focus
- No Cython/C compilation issues
- Focus on Python best practices
- Leverage scientific Python ecosystem (NumPy, SciPy)

### Dependency Awareness
- **benpy** provides multi-objective optimization
- **COBRApy** provides metabolic models
- Don't reimplement what dependencies provide
- Understand integration points

### Computational Biology Context
- Users are biologists and bioinformaticians
- Documentation should explain biological concepts
- Examples should be biologically motivated
- Results need biological interpretation

## Questions to Ask When Choosing Agents

1. **Is this a biological question?** → bioinformatics-agent
2. **Does this involve code implementation?** → dev_base
3. **Do we need tests?** → testing-agent
4. **Is documentation needed?** → docagent
5. **Are there security/quality concerns?** → security-agent

## Conclusion

Effective agent coordination leads to:
- Better biological validity
- Higher code quality
- Comprehensive testing
- Clear documentation
- Secure implementation

When in doubt, consult multiple agents and synthesize their recommendations, always prioritizing biological correctness for this computational biology package.
