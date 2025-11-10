# mocbapy Custom Agents

This directory contains custom GitHub Copilot agents specialized for different aspects of the mocbapy project.

## What are Custom Agents?

Custom agents are specialized AI assistants configured with specific expertise and instructions for particular domains. When you work on mocbapy, you can delegate tasks to these agents to leverage their specialized knowledge.

## Available Agents

### Core Development
- **[dev_base.md](dev_base.md)** - Python Developer
  - General Python development for mocbapy
  - COBRApy and benpy integration
  - New features and API design
  - Bug fixes and refactoring
  - Pure Python code (no Cython)

### Domain Expertise
- **[bioinformatics-agent.md](bioinformatics-agent.md)** - Bioinformatics & Metabolic Modeling Expert
  - Metabolic modeling and FBA expertise
  - Systems biology and microbial ecosystems
  - Model validation and interpretation
  - Biological result validation
  - COBRA model expertise
  - **UNIQUE TO MOCBAPY** - Not in benpy

### Quality & Testing
- **[testing-agent.md](testing-agent.md)** - Testing Expert
  - Unit and integration tests
  - Numerical validation
  - Biological model testing
  - COBRA model test fixtures
  - Test infrastructure

- **[security-agent.md](security-agent.md)** - Security & Code Quality Expert
  - Security audits and vulnerability fixes
  - Input validation
  - Code quality and best practices
  - Dependency security

### Documentation
- **[docagent.md](docagent.md)** - Documentation Expert
  - README and documentation updates
  - API documentation and examples
  - Installation guides
  - Biological concept explanations

## How to Use Agents

### Basic Usage

When working on an issue or task, consider which agent(s) are most appropriate:

1. **Identify the domain** of your task (development, biology, testing, docs, etc.)
2. **Consult the coordination guide** in [AGENT_COORDINATION.md](AGENT_COORDINATION.md)
3. **Use the appropriate agent** for the task
4. **Chain agents** for complex multi-domain tasks

### Examples

**Adding a new MO-FBA feature:**
→ Use `dev_base` → `testing-agent` → `docagent`

**Validating ecosystem model:**
→ Use `bioinformatics-agent`

**Improving documentation:**
→ Use `docagent`

**Adding test suite:**
→ Use `testing-agent`

**Security review:**
→ Use `security-agent`

**Biological interpretation of results:**
→ Use `bioinformatics-agent`

## Agent Coordination

For complex tasks involving multiple domains, see [AGENT_COORDINATION.md](AGENT_COORDINATION.md) for:
- Task-to-agent mapping
- Multi-agent workflows
- Best practices for agent coordination
- Common scenarios and examples

## Agent Design Principles

All agents in this directory follow these principles:

1. **Specialized Expertise**: Each agent has deep knowledge in a specific domain
2. **Clear Scope**: Well-defined responsibilities and boundaries
3. **Actionable Guidance**: Specific instructions and examples
4. **Project Context**: Understanding of mocbapy's architecture and goals
5. **Best Practices**: Following industry standards and project conventions
6. **Minimal Changes**: Making surgical, targeted modifications
7. **Biological Awareness**: Understanding the computational biology context

## Differences from benpy Agents

mocbapy agents are tailored to this project's unique characteristics:

### Not Included (compared to benpy)
- **No Cross-Platform Compiler Agent**: mocbapy is pure Python (no Cython/C compilation)
- **No CI/CD Agent**: (Can be added later if CI/CD workflows are created)
- **No Performance Agent**: (Can be added later if needed)

### Unique to mocbapy
- **Bioinformatics Agent**: Domain expert for metabolic modeling and systems biology
  - Understands FBA, MO-FBA, and COBRA models
  - Validates biological correctness
  - Interprets results in biological context
  - Essential for computational biology package

### Adapted from benpy
- **Python Developer**: Focuses on pure Python (no Cython), COBRApy/benpy integration
- **Testing Agent**: Emphasizes biological model validation and COBRA test fixtures
- **Documentation Agent**: Includes biological concept explanations
- **Security Agent**: Focuses on scientific Python security (input validation, numerical stability)

## Key Dependencies

mocbapy agents should understand these dependencies:
- **benpy**: Vector Linear Programming (VLP) solver - multi-objective optimization engine
- **COBRApy**: Constraint-Based Reconstruction and Analysis - metabolic model framework
- **optlang**: Solver abstraction layer (GLPK, Gurobi, CPLEX)
- **scipy**: Sparse matrix operations
- **numpy**: Numerical computations
- **pandas**: Data structures for FVA results
- **tqdm**: Progress bars

## File Naming Convention

- Agent files use descriptive names: `{domain}-agent.md`
- Main agents may use shorter names: `dev_base.md`, `docagent.md`
- All agents have YAML frontmatter with name and description

## YAML Frontmatter Format

```yaml
---
name: Agent Display Name
description: Brief description of agent's purpose
---
```

## Need Help?

- See [AGENT_COORDINATION.md](AGENT_COORDINATION.md) for guidance on choosing agents
- Read individual agent files for their specific capabilities
- Consult the main [copilot-instructions.md](../copilot-instructions.md) for project overview

## Contributing

When adding new agents:

1. Identify a clear, unmet need
2. Define the agent's scope and expertise
3. Write comprehensive instructions with examples
4. Update [AGENT_COORDINATION.md](AGENT_COORDINATION.md)
5. Update this README

Avoid creating agents with overlapping responsibilities or agents that are too narrow to be useful.

## Updating Agents

When updating agent configurations:

1. **Test changes**: Ensure agent instructions are clear and effective
2. **Maintain consistency**: Follow the established format and style
3. **Document expertise**: Clearly define what the agent can help with
4. **Provide examples**: Include concrete code examples where relevant
5. **Keep focused**: Don't make agents too broad or overlapping
6. **Consider biology**: Ensure biological context is preserved in updates
