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
        'ГAMMA', 'TETA'
    ])
    rs_name = tuple(['Init', 'PK', 'CA'])
    profiler = profiler.Profiling(
        input_file,
        geometrical_parameters,
        rs_name
    )
    profiler.get_data()
    profiler.coordinates()
    profiler.twist_law()
    profiler.save_curve_files()
