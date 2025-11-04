import math
import numpy as np
from scipy.integrate import solve_ivp
import sympy as sp
from sympy.parsing.sympy_parser import parse_expr
from sympy import lambdify
import re
import copy
import os
from d6py.default_initial_values import INITIAL_VALUES as DEFAULT_INITIAL_VALUES
from d6py.default_initial_values import INITIAL_ENERGY_SCALE as DEFAULT_INITIAL_ENERGY_SCALE

class D6Solver:
    def __init__(self, initial_energy=None, final_energy=None, beta_file=None, redundant_file=None):
        #_here = os.getcwd()
        #_here = os.path.join(_here, "d6py")
        _here = os.path.dirname(__file__)
        
        self.INITIAL_ENERGY_SCALE = initial_energy or DEFAULT_INITIAL_ENERGY_SCALE
        self.FINAL_ENERGY_SCALE = final_energy

        self.INITIAL_VALUES = copy.deepcopy(DEFAULT_INITIAL_VALUES)
        self.RUNNING_VARIABLES = {}
        self.Running_variables_symbols = {}
        self.ADDITIONAL_CONSTANTS = {
            "LoopParameter": 0
        }
        self.BETA_EXPRS = {}
        self.Beta_funcs_np = {}
        self.Key_to_index = {}

        self.Solution = None

        self.Beta_fns_file = beta_file or os.path.join(_here, "BetaSMEFT.txt")
        self.redundant_elements_file = redundant_file or os.path.join(_here, "Redundant_Elements.txt")

        self.redundant_lhs = []
        self.redundant_relations = []
        self.redundant_vars_rhs = []

        self.load_redundant_elements()

    def load_redundant_elements(self):
        with open(self.redundant_elements_file, "r") as f:
            r_lines = f.readlines()
        for line in r_lines:
            line = line.strip()
            if not line or line.startswith('#'):
                continue
            lhs, rhs = map(str.strip, line.split("=", 1))
            self.redundant_lhs.append(lhs)
            if rhs.startswith("conjugate(") and rhs.endswith(")"):
                inner_var = rhs[len("conjugate("):-1].strip()
                self.redundant_vars_rhs.append(inner_var)
                self.redundant_relations.append("conjugate")
            else:
                self.redundant_vars_rhs.append(rhs)
                self.redundant_relations.append("equal")

    def Set_value(self, variable_name, value):
        """
        Set the value of a variable in self.INITIAL_VALUES.
        Handles redundant variables without recursion.
        """
        variable_name = variable_name.replace(" ", "")  # remove spaces
        match = re.match(r"(\w+)(?:\[(.*?)\])?$", variable_name)
        if not match:
            raise ValueError(f"Invalid variable name: {variable_name}")
        key, indices = match.groups()
        if key not in self.INITIAL_VALUES:
            raise KeyError(f"Variable {key} not found in INITIAL_VALUES")
        
        # Convert indices to tuple if present
        idx = tuple(int(i) for i in indices.split(",")) if indices else None

        # Set primary variable
        if idx:
            self.INITIAL_VALUES[key][idx] = value
        else:
            self.INITIAL_VALUES[key] = value

        # --- LHS redundants ---
        if variable_name in self.redundant_lhs:
            redundant_index = self.redundant_lhs.index(variable_name)
            rhs_variable = self.redundant_vars_rhs[redundant_index]
            relation = self.redundant_relations[redundant_index]
            rhs_match = re.match(r"(\w+)(?:\[(.*?)\])?$", rhs_variable)
            rhs_key, rhs_indices = rhs_match.groups()
            rhs_idx = tuple(int(i) for i in rhs_indices.split(",")) if rhs_indices else None
            if relation == "equal":
                if rhs_idx:
                    self.INITIAL_VALUES[rhs_key][rhs_idx] = value
                else:
                    self.INITIAL_VALUES[rhs_key] = value
            elif relation == "conjugate":
                if rhs_idx:
                    self.INITIAL_VALUES[rhs_key][rhs_idx] = np.conjugate(value)
                else:
                    self.INITIAL_VALUES[rhs_key] = np.conjugate(value)

        # --- RHS redundants ---
        for i, rhs_var in enumerate(self.redundant_vars_rhs):
            if rhs_var == variable_name:
                lhs = self.redundant_lhs[i]
                relation = self.redundant_relations[i]
                lhs_match = re.match(r"(\w+)(?:\[(.*?)\])?$", lhs)
                lhs_key, lhs_indices = lhs_match.groups()
                lhs_idx = tuple(int(i) for i in lhs_indices.split(",")) if lhs_indices else None
                val_to_set = value if relation == "equal" else np.conjugate(value)
                if lhs_idx:
                    self.INITIAL_VALUES[lhs_key][lhs_idx] = val_to_set
                else:
                    self.INITIAL_VALUES[lhs_key] = val_to_set

    def Initialize_running_variables(self, running_variable_names, show_warnings=False):
        """
        Initialize self.RUNNING_VARIABLES from self.INITIAL_VALUES based on running_variable_names.
        running_variable_names: list of str, e.g., ['g', 'Gu[1,2]', 'Gu', etc.]
        show_warnings: bool, if True, print warnings for redundant variables found.
        """
        self.RUNNING_VARIABLES.clear()

        for expr in running_variable_names:
            expr = expr.replace(" ", "")  # remove spaces
            match = re.match(r"(\w+)(?:\[(.*?)\])?$", expr)
            if not match:
                raise ValueError(f"Invalid request: {expr}")
            key, indices = match.groups()
            if key not in self.INITIAL_VALUES:
                raise KeyError(f"Variable {key} not found in INITIAL_VALUES")
            
            value = self.INITIAL_VALUES[key]

            if indices:  # element access like Gu[1,2]
                if expr in self.redundant_lhs:
                    redundant_index = self.redundant_lhs.index(expr)
                    lhs = expr
                    rhs_var = self.redundant_vars_rhs[redundant_index]
                    relation = self.redundant_relations[redundant_index]
                    if show_warnings:
                        rhs = rhs_var if relation == "equal" else f"conjugate({rhs_var})"
                        print(f"Redundant variable found {lhs} = {rhs}. Using {rhs}.")
                    expr = rhs_var
                idx = tuple(int(i) for i in indices.split(","))
                value = value[idx]
                self.RUNNING_VARIABLES[expr] = value

            else:  # full variable
                if isinstance(value, np.ndarray):
                    for idx in np.ndindex(value.shape):
                        dict_key = f"{expr}[{','.join(map(str, idx))}]"
                        if dict_key not in self.redundant_lhs:
                            self.RUNNING_VARIABLES[dict_key] = value[idx]
                        else:
                            lhs = dict_key
                            redundant_index = self.redundant_lhs.index(dict_key)
                            rhs_var = self.redundant_vars_rhs[redundant_index]
                            relation = self.redundant_relations[redundant_index]
                            rhs = rhs_var if relation == "equal" else f"conjugate({rhs_var})"
                            if show_warnings:
                                print(f"Redundant variable found {lhs} = {rhs}. Using {rhs}.")
                else:
                    self.RUNNING_VARIABLES[expr] = value

    def load_beta_functions(self, Beta_fns_file=None):
        file_to_use = Beta_fns_file or self.Beta_fns_file
        self.BETA_EXPRS.clear()
        self.Running_variables_symbols.clear()

        for key, val in self.RUNNING_VARIABLES.items():
            key_mod = key.replace("[", "_").replace(",", "_").replace("]", "")
            self.Running_variables_symbols[key_mod] = sp.Symbol(key_mod)

        constant_variables = {}
        for key, val in self.INITIAL_VALUES.items():
            if isinstance(val, np.ndarray):
                for idx in np.ndindex(val.shape):
                    dict_key = f"{key}_{'_'.join(map(str, idx))}"
                    if dict_key not in self.Running_variables_symbols:
                        constant_variables[dict_key] = val[idx]
            else:
                if key not in self.Running_variables_symbols:
                    constant_variables[key] = val

        constant_variables.update(self.ADDITIONAL_CONSTANTS)
        constant_variables["Pi"] = np.pi
        constant_variables["Zeta_3"] = 1.2020569032
        constant_variables["Zeta_7"] = 1.036927755
        constant_variables['conjugate'] = sp.conjugate

        with open(file_to_use, "r") as f:
            lines = f.readlines()

        local_dict = {}
        local_dict.update(self.Running_variables_symbols)
        local_dict.update(constant_variables)

        temp_beta_exprs = {}
        for line in lines:
            line = line.strip()
            if not line or line.startswith('#'):
                continue
            if "=" in line:
                lhs, rhs = map(str.strip, line.split("=", 1))
                lhs_mod = lhs.replace("[", "_").replace(",", "_").replace("]", "")
                
                if lhs_mod.replace("Beta_", "") in self.Running_variables_symbols:
                    rhs_mod = rhs.replace("[", "_").replace(",", "_").replace("]", "")
                    expr = parse_expr(rhs_mod, local_dict=local_dict)
                    temp_beta_exprs[lhs] = expr

        ordered = {f"Beta_{k}": temp_beta_exprs[f"Beta_{k}"] for k in self.RUNNING_VARIABLES.keys()}
        self.BETA_EXPRS.update(ordered)
        self.update_numpy_functions()

    def update_numpy_functions(self):
        self.Beta_funcs_np.clear()
        for name, expr in self.BETA_EXPRS.items():
            self.Beta_funcs_np[name] = lambdify(tuple(self.Running_variables_symbols.values()), expr, "numpy")

    def ode_system(self, t, y): # (1/(16*pi^2)) * beta_y
        beta_y = np.zeros_like(y)
        for i in range(len(y)):
            beta_y[i] = list(self.Beta_funcs_np.values())[i](*y)
        return (1/(16*np.pi**2)) * beta_y

    def Run(self):
        y0 = []
        for key, val in self.RUNNING_VARIABLES.items():
            if isinstance(val, np.ndarray):
                y0.extend(val.flatten())
            else:
                y0.append(val)

        self.Solution = solve_ivp(
            self.ode_system,
            t_span=(math.log(self.INITIAL_ENERGY_SCALE), math.log(self.FINAL_ENERGY_SCALE)),
            y0=y0,
            method='RK45',
            rtol=1e-6,
            atol=1e-9,
            dense_output=True
        )

        self.Key_to_index.clear()
        i = 0
        for key, val in self.RUNNING_VARIABLES.items():
            if isinstance(val, np.ndarray):
                self.Key_to_index[key] = [i, i + val.size]
                i += val.size
            else:
                self.Key_to_index[key] = [i, i + 1]
                i += 1

    def Get_values(self, mu, variable_name, solution=None):
        if solution is None:
            solution = self.Solution
        t = np.log(mu)
        try:
            index_range = self.Key_to_index[variable_name]
            return solution.sol(t)[index_range[0]:index_range[1]].flatten()
        except KeyError:
            if variable_name in self.redundant_lhs:
                redundant_index = self.redundant_lhs.index(variable_name)
                lhs = variable_name
                rhs = self.redundant_vars_rhs[redundant_index]
                relation = self.redundant_relations[redundant_index]
                if relation == "equal":
                    return self.Get_values(mu, rhs, solution=solution)
                elif relation == "conjugate":
                    return self.Get_values(mu, rhs, solution=solution).conjugate()
            raise KeyError(f"Variable {variable_name} not found in RUNNING_VARIABLES.")
