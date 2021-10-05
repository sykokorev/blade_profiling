# -*- coding: utf-8 -*-

import pandas as pd
import itertools
import math

from _collections_abc import MutableMapping

import Parser as Pr


def _dict_unpack(d, parent_key, sep):
    for k, v, in d.items():
        new_key = str(parent_key) + sep + str(k) if parent_key else str(k)
        if isinstance(v, MutableMapping):
            yield from dict_unpack(v, new_key, sep=sep).items()
        else:
            yield new_key, v


def dict_unpack(d: MutableMapping, parent_key: str = '', sep: str = '.'):
    return dict(_dict_unpack(d, parent_key, sep))


def list_keys_generator(keys_tuple: tuple, sep: str = '.'):
    keys = []
    for key in keys_tuple:
        new_key = ''
        for el in key:
            sep_temp = (sep if key.index(el) < (len(key) - 1) else '')
            new_key += f"{el}" + sep_temp
        keys.append(new_key)
    return keys


def get_unique_keys(*args):
    unique_keys = set(args[0])
    for el in args:
        unique_keys.intersection_update(el)
    return tuple(unique_keys)


class Profiling:

    def __init__(self, in_file=None, parameters=None, rs_name=None):

        self.parameters = parameters
        self.rs_name = rs_name
        self.hub_coordinates = []
        self.shroud_coordinates = []
        self.stages_num = 0
        self.df = pd.read_csv(in_file, header=None, encoding="cp1251")
        self.data = {}

    def get_data(self):
        self.stages_num = Pr.number_of_stages(self.df)

        for stage in range(1, self.stages_num + 1):
            self.data[stage] = Pr.data_for_stage(self.df, stage)

    def _data_formation(self, sep='.'):

        flat_dict = dict_unpack(self.data, sep=sep)

        keys_parameters = itertools.product(
            list(range(1, self.stages_num + 1)),
            self.rs_name,
            self.parameters,
            repeat=1
        )

        keys_angles = itertools.product(
            list(range(1, self.stages_num + 1)),
            self.parameters
        )

        flat_dict_keys = tuple(flat_dict.keys())
        keys = tuple(list_keys_generator(tuple(keys_parameters)) + list_keys_generator(tuple(keys_angles)))
        unique_keys = get_unique_keys(keys, flat_dict_keys)

        data = {}

        for stage in list(range(1, self.stages_num + 1)):
            data[stage] = {}

        for key in unique_keys:
            key_high = int(key.split('.')[0])
            key_low = ''
            for k in key.split(sep)[1:]:
                if key.split(sep)[1:].index(k) < (len(key.split(sep)[1:]) - 1):
                    s = sep
                else:
                    s = ''
                key_low += k + s
            data[key_high][key_low] = flat_dict[key]

        return data

    def coordinates(self, dz_=1.2):
        data_frame = pd.DataFrame(self._data_formation())
        print(data_frame)
        z = 0
        dz = {}

        for stage in list(range(1, self.stages_num + 1)):
            dz[stage] = [
                dz_ * data_frame[stage]['CA.T'] * math.tan(math.radians(data_frame[stage]['AL1'])),
                dz_ * data_frame[stage]['PK.T'] * math.tan(math.radians(data_frame[stage]['BE2']))
            ]
            throat_gv = 0.5 * data_frame[stage]['CA.A/T'] * data_frame[stage]['CA.T'] * \
                        math.cos(math.radians(data_frame[stage]['AL1']))
            throat_rb = 0.5 * data_frame[stage]['PK.A/T'] * data_frame[stage]['PK.T'] * \
                        math.cos(math.radians(data_frame[stage]['BE2']))
            z_throat_gv = z + data_frame[stage]['CA.B'] - throat_gv
            z_throat_rb = z_throat_gv + throat_gv + dz[stage][1] + data_frame[stage]['PK.B'] - throat_rb
