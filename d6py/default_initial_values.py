import numpy as np

INITIAL_ENERGY_SCALE = 91.1876

INITIAL_VALUES = {
    # === SM parameters ===
    'g': 0.629787,
    'gp': 0.345367,
    'gs': 1.218232,
    'Lambda': 0.257736,
    'm2': 7812.5,
    'Theta': 0.0,
    'Theta_p': 0.0,
    'Theta_s': 0.0,

    # === Yukawas ===
    'Gu': np.array([
        [1.23231E-5, -1.64215E-3, 5.90635E-3],
        [2.84527E-6, 7.10724E-3, -4.18547E-2],
        [4.65426E-8, 3.08758E-4, 0.994858]
    ]),
    'Gd': np.array([
        [2.70195E-5, 0, 0],
        [0, 5.51888E-4, 0],
        [0, 0, 2.403012E-2]
    ]),
    'Ge': np.array([
        [2.93766E-6, 0, 0],
        [0, 6.07422E-4, 0],
        [0, 0, 1.02157E-2]
    ]),

    # === Bosonic operators ===
    'CG': 0.0,
    'CGtilde': 0.0,
    'CW': 0.0,
    'CWtilde': 0.0,
    'CH': 0.0,
    'CHbox': 0.0,
    'CHD': 0.0,
    'CHG': 0.0,
    'CHB': 0.0,
    'CHW': 0.0,
    'CHWB': 0.0,
    'CHGtilde': 0.0,
    'CHBtilde': 0.0,
    'CHWtilde': 0.0,
    'CHWtildeB': 0.0,

    # === Fermion-Higgs operators ===
    'CuH': np.zeros((3, 3)),
    'CdH': np.zeros((3, 3)),
    'CeH': np.zeros((3, 3)),
    'CeW': np.zeros((3, 3)),
    'CeB': np.zeros((3, 3)),
    'CuG': np.zeros((3, 3)),
    'CuW': np.zeros((3, 3)),
    'CuB': np.zeros((3, 3)),
    'CdG': np.zeros((3, 3)),
    'CdW': np.zeros((3, 3)),
    'CdB': np.zeros((3, 3)),

    # === Fermion-gauge-Higgs currents ===
    'CHl1': np.zeros((3, 3)),
    'CHl3': np.zeros((3, 3)),
    'CHe': np.zeros((3, 3)),
    'CHq1': np.zeros((3, 3)),
    'CHq3': np.zeros((3, 3)),
    'CHu': np.zeros((3, 3)),
    'CHd': np.zeros((3, 3)),
    'CHud': np.zeros((3, 3)),

    # === Four-fermion operators ===
    'Cll': np.zeros((3, 3, 3, 3)),
    'Cqq1': np.zeros((3, 3, 3, 3)),
    'Cqq3': np.zeros((3, 3, 3, 3)),
    'Clq1': np.zeros((3, 3, 3, 3)),
    'Clq3': np.zeros((3, 3, 3, 3)),
    'Cee': np.zeros((3, 3, 3, 3)),
    'Cuu': np.zeros((3, 3, 3, 3)),
    'Cdd': np.zeros((3, 3, 3, 3)),
    'Ceu': np.zeros((3, 3, 3, 3)),
    'Ced': np.zeros((3, 3, 3, 3)),
    'Cud1': np.zeros((3, 3, 3, 3)),
    'Cud8': np.zeros((3, 3, 3, 3)),
    'Cle': np.zeros((3, 3, 3, 3)),
    'Clu': np.zeros((3, 3, 3, 3)),
    'Cld': np.zeros((3, 3, 3, 3)),
    'Cqe': np.zeros((3, 3, 3, 3)),
    'Cqu1': np.zeros((3, 3, 3, 3)),
    'Cqu8': np.zeros((3, 3, 3, 3)),
    'Cqd1': np.zeros((3, 3, 3, 3)),
    'Cqd8': np.zeros((3, 3, 3, 3)),
    'Cledq': np.zeros((3, 3, 3, 3)),
    'Cquqd1': np.zeros((3, 3, 3, 3)),
    'Cquqd8': np.zeros((3, 3, 3, 3)),
    'Clequ1': np.zeros((3, 3, 3, 3)),
    'Clequ3': np.zeros((3, 3, 3, 3)),
    'Cduql': np.zeros((3, 3, 3, 3)),
    'Cqque': np.zeros((3, 3, 3, 3)),
    'Cqqql': np.zeros((3, 3, 3, 3)),
    'Cduue': np.zeros((3, 3, 3, 3)),

    # === Special ===
    'CllHH': np.zeros((3, 3)),
}