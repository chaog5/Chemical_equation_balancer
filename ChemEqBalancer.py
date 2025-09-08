import re
import sys
from collections import defaultdict, OrderedDict
from sympy import Matrix, lcm, nsimplify
from fractions import Fraction

# Simplified regex pattern
chem_eq_regex = re.compile(r'^(?P<reactants>.+?)\s*(?:->|→|=)\s*(?P<products>.+)$', re.IGNORECASE)

class ChemEqParser:
    def __init__(self):
        self.reactants = []
        self.products = []
        self.chem_species = []
        self.chem_species_counts = []
        self.all_elements = OrderedDict()
        self.parsed_chem_eq = OrderedDict()

    def validate_formula(self, formula):
        """Check for common input errors like using numbers where letters should be"""
        # Check if formula ends with a number (like H20 instead of H2O)
        if formula[-1].isdigit() and any(c.isalpha() for c in formula[:-1]):
            # Look for patterns like H20, O3H20, etc.
            elements = []
            i = 0
            while i < len(formula):
                if formula[i].isupper():
                    element_start = i
                    i += 1
                    while i < len(formula) and formula[i].islower():
                        i += 1
                    element = formula[element_start:i]
                    elements.append(element)
                else:
                    i += 1
            
            # If we found elements but the formula ends with a number, it might be a mistake
            if elements and formula[-1].isdigit():
                suggested = formula
                for common_error in [('0', 'O'), ('1', 'I')]:
                    if common_error[0] in formula:
                        suggested = suggested.replace(common_error[0], common_error[1])
                
                if suggested != formula:
                    raise ValueError(f"Invalid formula '{formula}'. Did you mean '{suggested}'?")

    def atom_count_in_formula(self, formula):
        # Validate formula first
        self.validate_formula(formula)
        
        if '·' in formula:
            parts = formula.split("·")
            counts = defaultdict(int)
            for part in parts:
                part_counts = self.atom_count_in_formula(part)
                for element, count in part_counts.items():
                    counts[element] += count
            return dict(counts)
        
        atoms = defaultdict(int)
        stack = []
        current_element = defaultdict(int)
        i = 0
        n = len(formula)

        while i < n:
            if formula[i] == "(":
                stack.append(current_element)
                current_element = defaultdict(int)
                i += 1
            elif formula[i] == ")":
                multiplier = 0
                i += 1
                while i < n and formula[i].isdigit():
                    multiplier = multiplier * 10 + int(formula[i])
                    i += 1
                multiplier = multiplier if multiplier else 1
                for element, count in current_element.items():
                    if stack:
                        stack[-1][element] += count * multiplier
                current_element = stack.pop() if stack else current_element
            else:
                if formula[i].isupper():
                    element_begins = i
                    i += 1
                    while i < n and formula[i].islower():
                        i += 1
                    element = formula[element_begins:i]
                    count_num = 0
                    while i < n and formula[i].isdigit():
                        count_num = count_num * 10 + int(formula[i])
                        i += 1
                    current_element[element] += count_num if count_num else 1
                else:
                    i += 1
        return dict(current_element)

    def equation_parser(self, equation):
        parts = re.split(r'\s*(?:->|→|=)\s*', equation.strip())
        if len(parts) != 2:
            raise ValueError('Invalid chemical equation. Use "->", "→", or "=" to separate reactants and products')
        
        self.reactants = [species.strip() for species in parts[0].split("+")]
        self.products = [species.strip() for species in parts[1].split("+")]
        self.chem_species = self.reactants + self.products

        self._get_atom_count()
        self._get_species_dicts()

        return self.parsed_chem_eq

    def _get_atom_count(self):
        self.chem_species_counts = []
        self.all_elements = OrderedDict()

        for species in self.chem_species:
            counts = self.atom_count_in_formula(species)
            self.chem_species_counts.append(counts)
            for element in counts.keys():
                if element not in self.all_elements:
                    self.all_elements[element] = None

    def _get_species_dicts(self):
        self.parsed_chem_eq = OrderedDict()
        num_reactants = len(self.reactants)

        for i, species in enumerate(self.chem_species):
            species_dict = OrderedDict()
            for element in self.all_elements:
                count = self.chem_species_counts[i].get(element, 0)
                if i >= num_reactants:
                    count *= -1
                species_dict[element] = count
            self.parsed_chem_eq[species] = species_dict

    def get_matrix(self):
        if not self.parsed_chem_eq:
            raise ValueError('Need to parse the chemical equation. Call equation_parser() first')
        
        elements = list(self.all_elements.keys())
        matrix = []

        for element in elements:
            row = []
            for species in self.chem_species:
                row.append(self.parsed_chem_eq[species].get(element, 0))
            matrix.append(row)

        return {
            "matrix": matrix,
            "elements": elements,
            "species": self.chem_species
        }

def get_stoic(equation, show_work=False):
    try:
        parser = ChemEqParser()
        parser.equation_parser(equation)
        matrix_data = parser.get_matrix()
        
        A = Matrix(matrix_data["matrix"])

        nullspace = A.nullspace()

        if not nullspace:
            if not show_work:
                print("Equation cannot be balanced. The nullspace is empty.")
                print("This usually means the equation is invalid or impossible to balance. If you have catalyst(s) in the equation, please remove them and try again.")
            
            if show_work:
                return None, matrix_data, None, None, None, None
            return None
    
        # Get the first nullspace vector
        vec = nullspace[0]
        
        # Convert to fractions for precise arithmetic
        simplified_vec = [nsimplify(elem) for elem in vec]

        # Find the least common multiple of denominators
        denominators = []
        for elem in simplified_vec:
            if hasattr(elem, 'as_numer_denom'):
                num, denom = elem.as_numer_denom()
                denominators.append(denom)
            else:
                denominators.append(1)
                
        multiplier = lcm(denominators) if denominators else 1
        
        raw_coeffs = [elem * multiplier for elem in simplified_vec]

        # Check if all coefficients are positive integers
        if all(coeff.is_integer and coeff > 0 for coeff in raw_coeffs):
            final_coeffs = [int(coeff) for coeff in raw_coeffs]
        else:
            if any(coeff < 0 for coeff in raw_coeffs):
                abs_coeffs = [abs(coeff) for coeff in raw_coeffs]
                if all(coeff.is_integer and coeff > 0 for coeff in abs_coeffs):
                    final_coeffs = [int(coeff) for coeff in abs_coeffs]
                else:
                    final_coeffs = None
            else:
                final_coeffs = None

        if show_work:
            return final_coeffs, matrix_data, vec, simplified_vec, multiplier, raw_coeffs
        else:
            return final_coeffs
    
    except Exception as e:
        if show_work:
            return None, None, None, None, None, None
        return None

def get_balanced_chem_eq(equation):
    try:
        coeffs = get_stoic(equation)
        if coeffs is None:
            return None
        
        parser = ChemEqParser()
        parser.equation_parser(equation)
        species = parser.chem_species

        reactants = []
        products = []
        arrow = "->"
        if "=" in equation:
            arrow = "="
        elif "→" in equation:
            arrow = "→"
        
        parts = equation.split(arrow)
        num_reactants = len(parts[0].split("+"))

        for i, (coeff, spec) in enumerate(zip(coeffs, species)):
            formatted_spec = f'{coeff if coeff != 1 else ""}{spec}'.strip()
            if i < num_reactants:
                reactants.append(formatted_spec)
            else:
                products.append(formatted_spec)
        
        balanced_eq = " + ".join(reactants) + f' {arrow} ' + " + ".join(products)

        return balanced_eq
    except ValueError as e:
        print(f"Input error: {e}")
        return None
    
def show_work(equation):
    coeffs, matrix_data, vec, simplified_vec, multiplier, raw_coeffs = get_stoic(equation, show_work=True)
    
    if matrix_data is None:
        print("Could not parse the equation.")
        return
    
    print(f"\nMatrix: {matrix_data['matrix']}")
    print(f"Elements: {matrix_data['elements']}")
    print(f"Species: {matrix_data['species']}")
    
    if vec is not None:
        print(f"Nullspace vector: {vec}")
        print(f"Simplified vector: {simplified_vec}")
        print(f"Multiplier: {multiplier}")
        print(f"Raw coefficients: {raw_coeffs}")
        print(f"Final coefficients: {coeffs}")
    else:
        print("No nullspace found - equation cannot be balanced")

def main():
    print('Enter an unbalanced chemical equation (e.g., "H2 + O2 -> H2O"). To quit, enter "q"\n')
    print('Note: Use letter "O" for oxygen, not the number "0"\n')

    last_balanced_eq = None

    while True:
        input_eq = input("\nEnter equation: ").strip()

        if input_eq.lower() == "q":
            print("\nGood-bye")
            sys.exit(0)
        
        if not input_eq:
            print('Please enter an unbalanced chemical equation or "q" to quit')
            continue
        
        if input_eq == 'show work':
            if last_balanced_eq:
                show_work(last_balanced_eq)
            else:
                print("No previous balanced equation to show work for.")
            continue

        last_balanced_eq = input_eq

        balanced_eq = get_balanced_chem_eq(input_eq)
        if balanced_eq:
            print(f"\nBalanced equation:\n{balanced_eq}")
            print(f'To show work, enter "show work". To quit, enter "q".')
        print("\n" + "="*50)
        



if __name__ == "__main__":
    main()