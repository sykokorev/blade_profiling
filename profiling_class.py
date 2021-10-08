# -*- coding: utf-8 -*-
import os

import pandas as pd
import itertools
import math


from _collections_abc import MutableMapping


import Parser as Pr


from common_class import CommonClass


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

        self.root_dir = os.getcwd()
        self.in_file = in_file
        self.parameters = parameters
        self.rs_name = rs_name
        self.hub = {}
        self.tip = {}
        self.z = {}
        self.r_tip = {}
        self.r_hub = {}
        self.beta = {}
        self.alfa = {}
        self.throat_width = {}
        self.stages_num = 0
        self.data_frame = pd.read_csv(in_file, header=None, encoding="cp1251")
        self.data = {}

    def get_data(self):
        self.stages_num = Pr.number_of_stages(self.data_frame)

        for stage in range(1, self.stages_num + 1):
            self.data[stage] = Pr.data_for_stage(self.data_frame, stage)

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
        df = pd.DataFrame(self._data_formation())

        d_throat_tip = {}
        d_throat_hub = {}
        z0 = 0
        dz = {}
        z_coordinates = {}
        d_tip = {}
        d_hub = {}

        for stage in list(range(1, self.stages_num + 1)):

            d_tip[stage] = []
            d_hub[stage] = []
            d_throat_tip[stage] = []
            d_throat_hub[stage] = []

            self.z[stage] = []

            dz[stage] = [
                dz_ * df[stage]['CA.T'] *
                math.tan(math.radians(df[stage]['AL1'])),
                dz_ * df[stage]['PK.T'] *
                math.tan(math.radians(df[stage]['BE2']))
            ]

            if stage > 1:
                z0 += df[stage]['CA.L'] * math.sin(math.radians(df[stage]['CA.ГAMMA']))\
                      + dz[stage][0] + df[stage]['PK.L'] * \
                      math.sin(math.radians(df[stage]['PK.ГAMMA'])) + dz[stage][1]

            throat_gv = 0.5 * df[stage]['CA.A/T'] * df[stage]['CA.T'] * \
                        math.cos(math.radians(df[stage]['AL1']))
            throat_rb = 0.5 * df[stage]['PK.A/T'] * df[stage]['PK.T'] * \
                        math.cos(math.radians(df[stage]['BE2']))

            z_throat_gv = z0 + df[stage]['CA.L'] * \
                          math.sin(math.radians(df[stage]['CA.ГAMMA'])) - throat_gv
            z1 = z0 + df[stage]['CA.L'] * math.sin(math.radians(df[stage]['CA.ГAMMA']))
            z2 = z1 + dz[stage][0]
            z_throat_rb = z2 + df[stage]['PK.L'] * \
                          math.sin(math.radians(df[stage]['PK.ГAMMA'])) - throat_rb
            z3 = z2 + df[stage]['PK.L'] * math.sin(math.radians(df[stage]['PK.ГAMMA']))

            z_coordinates[stage] = [z0, z_throat_gv, z1, z2, z_throat_rb, z3]

            # Line equation
            d_throat_tip[stage].append(df[stage]['Init.D1'] + df[stage]['CA.HЛ'])
            d_throat_tip[stage].append(df[stage]['Init.D2'] + df[stage]['PK.HЛ'])
            d_throat_hub[stage].append(df[stage]['Init.D1'] - df[stage]['CA.HЛ'])
            d_throat_hub[stage].append(df[stage]['Init.D2'] - df[stage]['PK.HЛ'])

            self.z[stage] = tuple([z0, z1, z2, z3])

        for stage in z_coordinates.keys():
            self.r_hub[stage] = []
            self.r_tip[stage] = []

            for z in z_coordinates[stage]:

                temp = (z - z_coordinates[stage][1])
                temp /= (z_coordinates[stage][4] - z_coordinates[stage][1])
                temp *= (d_throat_tip[stage][1] - d_throat_tip[stage][0])
                temp += d_throat_tip[stage][0]
                d_tip[stage].append(temp)

                if z in self.z[stage]:
                    self.r_tip[stage].append(temp/2)
                temp = (z - z_coordinates[stage][1])
                temp /= (z_coordinates[stage][4] - z_coordinates[stage][1])
                temp *= (d_throat_hub[stage][1] - d_throat_hub[stage][0])
                temp += d_throat_hub[stage][0]
                d_hub[stage].append(temp)

                if z in self.z[stage]:
                    self.r_hub[stage].append(temp/2)

            self.tip[stage] = zip(z_coordinates[stage], d_tip[stage])
            self.hub[stage] = zip(z_coordinates[stage], d_hub[stage])

    def twist_law(self):
        df = pd.DataFrame(self._data_formation())
        print(df)

        for stage in list(range(1, self.stages_num+1)):

            # Alfa GV
            # alfa_0/1[0] - hub, alfa_0/1[1] - mid, alfa_0/1[2] - shroud
            alfa_0_mid = math.radians(df[stage]['AL0'])
            alfa_1_mid = math.radians(df[stage]['AL1'])
            r_mid_gv = (self.r_tip[stage][1] + self.r_hub[stage][1]) / 2
            alfa_1 = tuple([
                math.degrees(alfa_1_mid * (self.r_hub[stage][1] / r_mid_gv)),
                math.degrees(alfa_1_mid),
                math.degrees(alfa_1_mid * (self.r_tip[stage][1] / r_mid_gv)),
            ])
            alfa_0 = tuple([
                math.degrees(alfa_0_mid * (self.r_hub[stage][1] / r_mid_gv)),
                math.degrees(alfa_0_mid),
                math.degrees(alfa_0_mid * (self.r_tip[stage][1] / r_mid_gv)),
            ])
            self.alfa[stage] = [alfa_0, alfa_1]

            # Beta RB
            # beta_1/2[0] - hub, beta_1/2[1] - mid, beta_1/2[2] - shroud
            beta_1 = tuple([
                df[stage]['BE1K'],
                df[stage]['BE1'],
                df[stage]['BE1П']
            ])
            beta_2_mid = math.radians(df[stage]['BE2'])
            r_mid_rb = (self.r_tip[stage][3] + self.r_hub[stage][3]) / 2
            beta_2 = tuple([
                math.degrees(math.tan(beta_2_mid) * (r_mid_rb / self.r_hub[stage][3])),
                math.degrees(beta_2_mid),
                math.degrees(math.tan(beta_2_mid) * (r_mid_rb / self.r_tip[stage][3]))
            ])
            self.beta[stage] = [beta_1, beta_2]

            # Throat width
            # throat_width[stage][GV(hub, mid, tip), RB(hub, mid, tip)]
            self.throat_width[stage] = [tuple([
                math.sin(math.radians(alfa_1[0])) * 2 *
                math.pi * self.r_hub[stage][1] / df[stage]['CA.ZЛ'],
                math.sin(math.radians(alfa_1[1])) * 2 *
                math.pi * r_mid_gv / df[stage]['CA.ZЛ'],
                math.sin(math.radians(alfa_1[2])) * 2 *
                math.pi * self.r_tip[stage][1] / df[stage]['CA.ZЛ']
            ]), tuple([
                math.sin(math.radians(beta_2[0])) * 2 *
                math.pi * self.r_hub[stage][3] / df[stage]['PK.ZЛ'],
                math.sin(math.radians(beta_2[1])) * 2 *
                math.pi * r_mid_rb / df[stage]['PK.ZЛ'],
                math.sin(math.radians(beta_2[2])) * 2 *
                math.pi * self.r_tip[stage][3] / df[stage]['PK.ZЛ']
            ])]

    def ab_data(self):
        pass

    def save_curve_files(self):
        dir_handle = CommonClass()
        curves_dir = os.path.join(self.root_dir, 'curves')
        log_msg, is_created = dir_handle.create_dir(curves_dir)
        if is_created:
            out_file_name_shroud = os.path.splitext(os.path.split(self.in_file)[1])[0] + '_shroud.dat'
            out_file_name_hub = os.path.splitext(os.path.split(self.in_file)[1])[0] + '_hub.dat'
            out_file_name_shroud = os.path.join(curves_dir, out_file_name_shroud)
            out_file_name_hub = os.path.join(curves_dir, out_file_name_hub)
