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

    # AutoBlade parameters
    twist_law = 1

    throat_base_additional_args = {
        'S1_SPAN': 0,
        'S2_SPAN': 0.5,
        'S3_SPAN': 1,
        # 'S1_SS_LE_HALF_WEDGE_ANGLE': 4.5,
        # 'S2_SS_LE_HALF_WEDGE_ANGLE': 4.5,
        # 'S3_SS_LE_HALF_WEDGE_ANGLE': 4.5,
        # 'S1_PS_LE_HALF_WEDGE_ANGLE': 4.5,
        # 'S2_PS_LE_HALF_WEDGE_ANGLE': 4.5,
        # 'S3_PS_LE_HALF_WEDGE_ANGLE': 4.5,
        # 'S1_DEV_ANGLE': 11,
        # 'S2_DEV_ANGLE': 11,
        # 'S3_DEV_ANGLE': 11,
        # 'S1_TE_WEDGE_ANGLE': 6.5,
        # 'S2_TE_WEDGE_ANGLE': 6.5,
        # 'S3_TE_WEDGE_ANGLE': 6.5,
        'S1_TPS_1': 0.06,
        'S1_TPS_2': 0.05,
        'S1_TPS_3': 0.01,
        'S1_TPS_4': -0.05,
        'S1_TSS_1': 0.07,
        'S1_TSS_2': 0.1,
        'S1_TSS_3': 0.18,
        'S1_TSS_4': 0.35,
        'S2_TPS_1': 0.06,
        'S2_TPS_2': 0.05,
        'S2_TPS_3': 0.01,
        'S2_TPS_4': -0.05,
        'S2_TSS_1': 0.07,
        'S2_TSS_2': 0.1,
        'S2_TSS_3': 0.18,
        'S2_TSS_4': 0.35,
        'S3_TPS_1': 0.06,
        'S3_TPS_2': 0.05,
        'S3_TPS_3': 0.01,
        'S3_TPS_4': -0.05,
        'S3_TSS_1': 0.07,
        'S3_TSS_2': 0.1,
        'S3_TSS_3': 0.18,
        'S3_TSS_4': 0.35,
        'LEAN_BETA': 0,
        'BLADE_ROTATION': 0
    }

    rs_name = tuple(['Init', 'PK', 'CA'])
    profiler = profiler.Profiling(
        input_file,
        geometrical_parameters,
        rs_name
    )
    profiler.get_data()
    profiler.save_curve_files()
    profiler.generate_param_file(**throat_base_additional_args)
