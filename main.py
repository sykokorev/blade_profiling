# -*- coding: utf-8 -*-

import os

import profiling_class as profiler

if __name__ == "__main__":
    root_dir = os.getcwd()
    data_dir = os.path.join(root_dir, "data")

    input_file = os.path.join(data_dir, "St_base_m90.REZ")
    geometrical_parameters = tuple([
        'AL0', 'AL1', 'AL2', 'BE2', 'BE1K', 'BE1П', 'BE1',
        'HЛ', 'T', 'A/T', 'D1', 'D2', 'B', 'A/T', 'ZЛ', 'L',
        'ГAMMA', 'TETA', 'HAД', 'PO', 'POK', 'FИ', 'V1', 'LU',
        'AL2'
    ])

    # Laws
    laws = {
        'stacking_law': 0,      # 0 - Center of gravity; 1 - LE, 2 - TE
        'surface_setup': 1,     # 0 - Conical, 1 - Cylindrical
        'twist_law': 1,         # 0 - Constant circulation, 1 - Constant alfa1
        'profile_law': 0        # 0 - Throat base, 1 - Legacy
    }

    additional_args = {
        'S1_SPAN': 0,
        'S2_SPAN': 0.5,
        'S3_SPAN': 1,
        'LEAN_BETA': 0,
        'BLADE_ROTATION': 0
    }

    rs_name = tuple(['Init', 'PK', 'CA'])
    profiler = profiler.Profiling(
        input_file,
        geometrical_parameters,
        rs_name,
        **laws
    )
    profiler.get_data()
    profiler.save_curve_files()
    profiler.generate_param_file(**additional_args)
