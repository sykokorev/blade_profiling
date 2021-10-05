# -*- coding: utf-8 -*-
import profiling_class as profiler


if __name__ == "__main__":
    input_file = r"./data/St_base_m90.REZ"
    geometrical_parameters = tuple([
        'AL0', 'AL1', 'AL2', 'BE2', 'BE1K', 'BE1П', 'BE1',
        'HЛ', 'T', 'A/T', 'D1', 'D2', 'B', 'A/T'
    ])
    rs_name = tuple(['Init', 'PK', 'CA'])
    profiler = profiler.Profiling(
        input_file,
        geometrical_parameters,
        rs_name
    )
    profiler.get_data()
    profiler.coordinates()
    # profiler._data_formation()
