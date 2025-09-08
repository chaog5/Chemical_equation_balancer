# Chemical_equation_balancer
A Python script that help chemists to balance chemical equations using linear algebra and matrix operations.

# Features:
  1. Automatic Balancing: Uses matrix nullspace calculations to determine stoichiometric coefficients
  2. Error Detection: Identifies common input errors (e.g., using numbers instead of letters)
  3. Work Display: Option to show the mathematical steps behind the balancing process
  4. Multiple Arrow Support: Works with ->, →, or = as equation separators
  5. Hydrate Support: Handles hydrated compounds (e.g., CuSO4·5H2O)

# Requirements:
  1. Python 3.x
  2. sympy

```bash

pip install sympy

```

# Usage

## Run the script:

```bash

python ChemEqBalancer.py

```

## Enter chemical equations in the format:

```text

Reactant1 + Reactant2 + ... -> Product1 + Product2 + ...

```

# Examples:
Input:

```text

H2 + O2 -> H2O

```

Output:

```text

Balanced equation:
2H2 + O2 -> 2H2O

```

Input:

```text

Fe + O2 = Fe2O3

```

Output:

```text

Balanced equation:
2H2 + O2 = 2H2O

```

Input:

```text
Al + H2SO4 → Al2(SO4)3 + H2

```

Output:

```text

Balanced equation:
2Al + 3H2SO4 → Al2(SO4)3 + 3H2

```

Input:

```text
CuSO4·5H2O -> CuSO4 + H2O

```

Output:

```text
Balanced equation:
CuSO4·5H2O -> CuSO4 + 5H2O

```

# Special Commands

- `q` - Quit the program
- `show work` - Display the mathematical steps for the last balanced equation

# How It Works:

  1. Parsing: The equation is split into reactants and products
  2. Element Counting: Each chemical formula is broken down into its constituent elements
  3. Matrix Formation: Creates a matrix where rows represent elements and columns represent compounds
  4. Nullspace Calculation: Uses linear algebra to find the solution to the homogeneous system
  5. Coefficient Adjustment: Converts fractional coefficients to whole numbers

# Supported Notation:

  1. Element symbols: H, He, Li, Be, etc. Essentially all standard chemical symbles.
  2. Subscripts: H₂O, CO₂, CH₄
  3. Parentheses: (NH₄)₂SO₄, Ca(OH)₂
  4. Hydrates: CuSO₄·5H₂O
  5. Multiple arrow types: ->, →, =

# Limitations:

Does not handle redox reactions specially
Cannot balance equations with catalysts. However, you can simply remove catalysts from your chemical equation, and the rest of reactants and products will be balanced for you.
May struggle with very complex equations

# Error Messages:

## The script provides helpful error messages for:

Using numbers instead of letters (e.g., "H20" instead of "H₂O")
Invalid chemical formulas
Equations that cannot be balanced

# Dependencies:

sympy: For matrix operations and nullspace calculation

# License:

This project is open source and available under the MIT License.

# Author:
Chao Gao, PhD (in chemistry)
[![GitHub](https://img.shields.io/badge/GitHub-@chaog5-blue?style=flat&logo=github)](https://github.com/chaog5)
