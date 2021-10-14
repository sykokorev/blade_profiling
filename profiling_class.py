# -*- coding: utf-8 -*-
import os

import pandas as pd
import itertools
import math

from _collections_abc import MutableMapping

import Parser as Pr

from common_class import CommonClass

pd.set_option('display.max_rows', None)
pd.set_option('display.max_columns', None)
pd.set_option('display.max_colwidth', None)


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

        self.dir_handle = CommonClass()

        self.root_dir = os.getcwd()
        self.names = ['param', 'rs']
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
        self.ab_data = pd.DataFrame()

    @staticmethod
    def line_equation(x1, y1, x2, y2, yi):
        xi = x2 - ((y2 - yi) / (y2 - y1)) * (x2 - x1)
        return xi

    @staticmethod
    def slope_angle(x1, y1, x2, y2):
        slope = math.degrees(math.atan((y2 - y1) / (x2 - x1)))
        return slope

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

    def _coordinates(self, dz_=1.2, sweep_beta=0, stacking_law=0):
        # Stacking law 0 - Center of gravity
        # 1 - Leading edge
        # 2 - Trailing edge
        # in any other case - Trailing edge
        df = pd.DataFrame(self._data_formation())
        d_throat_tip = {}
        d_throat_hub = {}
        z0 = 0
        dz = {}
        z_coordinates = {}
        b1 = {}

        zth = {}
        r_tip = {}
        r_hub = {}

        d_tip = {}
        d_hub = {}
        ref_trace = {}
        sweep = {}

        for stage in list(range(1, self.stages_num + 1)):

            d_tip[stage] = []
            d_hub[stage] = []
            d_throat_tip[stage] = []
            d_throat_hub[stage] = []

            self.z[stage] = []

            dz[stage] = [
                dz_ * df[stage]['CA.T'] * 1000 *
                math.tan(math.radians(df[stage]['AL1'])),
                dz_ * df[stage]['PK.T'] * 1000 *
                math.tan(math.radians(df[stage]['BE2']))
            ]

            if stage > 1:
                z0 += (df[stage]['CA.L'] * 1000 * math.sin(math.radians(df[stage]['CA.ГAMMA']))
                      + dz[stage][0] + df[stage]['PK.L'] * 1000 *
                      math.sin(math.radians(df[stage]['PK.ГAMMA'])) + dz[stage][1])

            throat_gv = (0.5 * df[stage]['CA.A/T'] * df[stage]['CA.T'] * 1000 *
                        math.cos(math.radians(df[stage]['AL1'])))
            throat_rb = 0.5 * df[stage]['PK.A/T'] * df[stage]['PK.T'] * \
                        math.cos(math.radians(df[stage]['BE2']))

            z_throat_gv = z0 + df[stage]['CA.L'] * 1000 * \
                          math.sin(math.radians(df[stage]['CA.ГAMMA'])) - throat_gv
            z1 = z0 + df[stage]['CA.L'] * 1000 * math.sin(math.radians(df[stage]['CA.ГAMMA']))
            z2 = z1 + dz[stage][0]
            z_throat_rb = z2 + df[stage]['PK.L'] * 1000 * \
                          math.sin(math.radians(df[stage]['PK.ГAMMA'])) - throat_rb
            z3 = z2 + df[stage]['PK.L'] * 1000 * math.sin(math.radians(df[stage]['PK.ГAMMA']))

            z_coordinates[stage] = [z0, z_throat_gv, z1, z2, z_throat_rb, z3]

            d_throat_tip[stage].append((df[stage]['Init.D1'] + df[stage]['CA.HЛ']) * 1000)
            d_throat_tip[stage].append((df[stage]['Init.D2'] + df[stage]['PK.HЛ']) * 1000)
            d_throat_hub[stage].append((df[stage]['Init.D1'] - df[stage]['CA.HЛ']) * 1000)
            d_throat_hub[stage].append((df[stage]['Init.D2'] - df[stage]['PK.HЛ']) * 1000)

            self.z[stage] = tuple([z0, z1, z2, z3])
            zth[stage] = tuple([z0-z0*0.1, z1+z1*0.1, z2-z2*0.1, z3+z3*0.1])
            r_tip[stage] = tuple([
                self.line_equation(
                    d_throat_tip[stage][0] / 2,
                    z_coordinates[stage][1],
                    d_throat_tip[stage][1] / 2,
                    z_coordinates[stage][4],
                    zth[stage][0]
                ),
                self.line_equation(
                    d_throat_tip[stage][0] / 2,
                    z_coordinates[stage][1],
                    d_throat_tip[stage][1] / 2,
                    z_coordinates[stage][4],
                    zth[stage][1]
                ),
                self.line_equation(
                    d_throat_tip[stage][0] / 2,
                    z_coordinates[stage][1],
                    d_throat_tip[stage][1] / 2,
                    z_coordinates[stage][4],
                    zth[stage][2]
                ),
                self.line_equation(
                    d_throat_tip[stage][0] / 2,
                    z_coordinates[stage][1],
                    d_throat_tip[stage][1] / 2,
                    z_coordinates[stage][4],
                    zth[stage][3]
                )
            ])
            r_hub[stage] = tuple([
                self.line_equation(
                    d_throat_hub[stage][0] / 2,
                    z_coordinates[stage][1],
                    d_throat_hub[stage][1] / 2,
                    z_coordinates[stage][4],
                    zth[stage][0]
                ),
                self.line_equation(
                    d_throat_hub[stage][0] / 2,
                    z_coordinates[stage][1],
                    d_throat_hub[stage][1] / 2,
                    z_coordinates[stage][4],
                    zth[stage][1]
                ),
                self.line_equation(
                    d_throat_hub[stage][0] / 2,
                    z_coordinates[stage][1],
                    d_throat_hub[stage][1] / 2,
                    z_coordinates[stage][4],
                    zth[stage][2]
                ),
                self.line_equation(
                    d_throat_hub[stage][0] / 2,
                    z_coordinates[stage][1],
                    d_throat_hub[stage][1] / 2,
                    z_coordinates[stage][4],
                    zth[stage][3]
                )
            ])

        for stage in z_coordinates.keys():
            self.r_hub[stage] = []
            self.r_tip[stage] = []

            for z in z_coordinates[stage]:

                temp = self.line_equation(
                    d_throat_tip[stage][0] / 2,
                    z_coordinates[stage][1],
                    d_throat_tip[stage][1] / 2,
                    z_coordinates[stage][4],
                    z
                )

                d_tip[stage].append(temp * 2)

                if z in self.z[stage]:
                    self.r_tip[stage].append(temp)

                temp = self.line_equation(
                    d_throat_hub[stage][0] / 2,
                    z_coordinates[stage][1],
                    d_throat_hub[stage][1] / 2,
                    z_coordinates[stage][4],
                    z
                )
                d_hub[stage].append(temp * 2)

                if z in self.z[stage]:
                    self.r_hub[stage].append(temp)

            self.tip[stage] = zip(z_coordinates[stage], d_tip[stage])
            self.hub[stage] = zip(z_coordinates[stage], d_hub[stage])
            ref_trace[stage] = tuple([
                self.z[stage][0], self.z[stage][0],
                self.r_tip[stage][0],
                self.r_hub[stage][0],
                self.slope_angle(self.z[stage][0], self.r_tip[stage][0],
                                 self.z[stage][1], self.r_tip[stage][1]),
                self.slope_angle(self.z[stage][0], self.r_hub[stage][0],
                                 self.z[stage][1], self.r_hub[stage][1]),
                self.z[stage][2], self.z[stage][2],
                self.r_tip[stage][2],
                self.r_hub[stage][2],
                self.slope_angle(self.z[stage][2], self.r_tip[stage][2],
                                 self.z[stage][3], self.r_tip[stage][3]),
                self.slope_angle(self.z[stage][2], self.r_hub[stage][2],
                                 self.z[stage][3], self.r_hub[stage][3]),
            ])
            if stacking_law == 0:
                sweep[stage] = [self.z[stage][0] + (self.z[stage][1] - self.z[stage][0]) / 2,
                                sweep_beta,
                                self.z[stage][2] + (self.z[stage][3] - self.z[stage][2]) / 2,
                                sweep_beta]
            elif stacking_law == 1:
                sweep[stage] = [self.z[stage][0],
                                sweep_beta,
                                self.z[stage][2],
                                sweep_beta]
            elif stacking_law == 2:
                sweep[stage] = [self.z[stage][1],
                                sweep_beta,
                                self.z[stage][3],
                                sweep_beta]
            else:
                sweep[stage] = [self.z[stage][1],
                                sweep_beta,
                                self.z[stage][3],
                                sweep_beta]

            b1[stage] = (
                self.z[stage][1] - self.z[stage][0],
                self.z[stage][1] - self.z[stage][0],
                self.z[stage][1] - self.z[stage][0],
                self.z[stage][3] - self.z[stage][2],
                self.z[stage][3] - self.z[stage][2],
                self.z[stage][3] - self.z[stage][2]
            )

        arrays = [
            ['HUB_R1', 'HUB_R2']*2,
            ['gv', 'gv', 'rb', 'rb']
        ]
        index = pd.MultiIndex.from_tuples(
            list(zip(*arrays)), names=self.names
        )
        coordinates = pd.DataFrame(r_hub, index=index)

        arrays = [
            ['SHROUD_R1', 'SHROUD_R2']*2,
            ['gv', 'gv', 'rb', 'rb']
        ]
        index = pd.MultiIndex.from_tuples(
            list(zip(*arrays)), names=self.names
        )
        coordinates = coordinates.append(pd.DataFrame(r_tip, index=index))

        arrays = [
            ['HUB_Z1', 'HUB_Z2']*2,
            ['gv', 'gv', 'rb', 'rb']
        ]
        index = pd.MultiIndex.from_tuples(
            list(zip(*arrays)), names=self.names
        )
        coordinates = coordinates.append(pd.DataFrame(zth, index=index))

        arrays = [
            ['SHROUD_Z1', 'SHROUD_Z2']*2,
            ['gv', 'gv', 'rb', 'rb']
        ]
        index = pd.MultiIndex.from_tuples(
            list(zip(*arrays)), names=self.names
        )
        coordinates = coordinates.append(pd.DataFrame(zth, index=index))

        arrays = [
            ['REF_TRACE_TIP_Z', 'REF_TRACE_HUB_Z',
             'REF_TRACE_TIP_R', 'REF_TRACE_HUB_R',
             'REF_TRACE_TIP_ALPHA', 'REF_TRACE_HUB_ALPHA'] * 2,
            ['gv']*6 + ['rb']*6
        ]
        index = pd.MultiIndex.from_tuples(
            list(zip(*arrays)), names=self.names
        )
        coordinates = coordinates.append(pd.DataFrame(ref_trace, index=index))

        arrays = [
            ['Z_STACKING_SWEEP', 'SWEEP_BETA']*2,
            ['gv']*2 + ['rb']*2
        ]
        index = pd.MultiIndex.from_tuples(
            list(zip(*arrays)), names=self.names
        )
        coordinates = coordinates.append(pd.DataFrame(sweep, index=index))

        arrays = [
            ['S1_DZ_PRIM', 'S2_DZ_PRIM', 'S3_DZ_PRIM']*2,
            ['gv']*3 + ['rb']*3
        ]

        index = pd.MultiIndex.from_tuples(
            list(zip(*arrays)), names=self.names
        )

        coordinates = coordinates.append(pd.DataFrame(b1, index=index))

        return coordinates

    def _twist_law(self, law=1):
        # law: 1 - twisted by constant circulation law;
        # law: 2 - twisted by constant reactivity
        df = pd.DataFrame(self._data_formation())

        if law == 1:
            for stage in list(range(1, self.stages_num + 1)):
                # Alfa GV
                # alfa_0/1[0] - hub, alfa_0/1[1] - mid, alfa_0/1[2] - shroud
                alfa_0_mid = math.radians(df[stage]['AL0'])
                alfa_1_mid = math.radians(df[stage]['AL1'])
                r_mid_gv = (self.r_tip[stage][1] + self.r_hub[stage][1]) / 2

                alfa_1 = [
                    math.degrees(alfa_1_mid * (self.r_hub[stage][1] / r_mid_gv)),
                    math.degrees(alfa_1_mid),
                    math.degrees(alfa_1_mid * (self.r_tip[stage][1] / r_mid_gv))
                ]

                if stage == 1:
                    alfa_0 = [float(90), float(90), float(90)]
                else:
                    alfa_0 = [
                        math.degrees(alfa_0_mid * (self.r_hub[stage][1] / r_mid_gv)),
                        math.degrees(alfa_0_mid),
                        math.degrees(alfa_0_mid * (self.r_tip[stage][1] / r_mid_gv)),
                    ]

                self.alfa[stage] = alfa_0 + alfa_1

                # Beta RB
                # beta_1/2[0] - hub, beta_1/2[1] - mid, beta_1/2[2] - shroud
                beta_1 = [
                    df[stage]['BE1K'],
                    df[stage]['BE1'],
                    df[stage]['BE1П']
                ]
                beta_2_mid = math.radians(df[stage]['BE2'])
                r_mid_rb = (self.r_tip[stage][3] + self.r_hub[stage][3]) / 2
                beta_2 = [
                    math.degrees(math.tan(beta_2_mid) * (r_mid_rb / self.r_hub[stage][3])),
                    math.degrees(beta_2_mid),
                    math.degrees(math.tan(beta_2_mid) * (r_mid_rb / self.r_tip[stage][3]))
                ]
                self.beta[stage] = beta_1 + beta_2

                # Throat width
                # throat_width[stage][GV(hub, mid, tip), RB(hub, mid, tip)]
                self.throat_width[stage] = [
                    math.sin(math.radians(alfa_1[0])) * 2 *
                    math.pi * self.r_hub[stage][1] / df[stage]['CA.ZЛ'],

                    # math.sin(math.radians(alfa_1[1])) * 2 *
                    # math.pi * r_mid_gv / df[stage]['CA.ZЛ'],
                    math.sin(math.radians(alfa_1[1])) * 2 *
                    math.pi * r_mid_gv / df[stage]['CA.ZЛ'],

                    math.sin(math.radians(alfa_1[2])) * 2 *
                    math.pi * self.r_tip[stage][1] / df[stage]['CA.ZЛ'],

                    math.sin(math.radians(beta_2[0])) * 2 *
                    math.pi * self.r_hub[stage][3] / df[stage]['PK.ZЛ'],

                    math.sin(math.radians(beta_2[1])) * 2 *
                    math.pi * r_mid_rb / df[stage]['PK.ZЛ'],

                    math.sin(math.radians(beta_2[2])) * 2 *
                    math.pi * self.r_tip[stage][3] / df[stage]['PK.ZЛ']
                ]

        elif law == 2:
            pass

        alfa_ab = {}
        for stage in self.alfa.keys():
            alfa_ab[stage] = [(x - 90) for x in self.alfa[stage][:3]]
            alfa_ab[stage] += [(90 - x) for x in self.alfa[stage][3:]]

        beta_ab = {}
        for stage in self.beta.keys():
            beta_ab[stage] = [(90 - x) for x in self.beta[stage][:3]]
            beta_ab[stage] += [(x - 90) for x in self.beta[stage][3:]]

        arrays = [
            ['S1_CAMBER_BETA1', 'S2_CAMBER_BETA1', 'S3_CAMBER_BETA1',
             'S1_CAMBER_BETA2', 'S2_CAMBER_BETA2', 'S3_CAMBER_BETA2'],
            ['gv'] * 6
        ]
        index = pd.MultiIndex.from_tuples(
            list(zip(*arrays)), names=self.names
        )
        angles_throat = pd.DataFrame(alfa_ab, index=index)

        arrays = [
            ['S1_CAMBER_BETA1', 'S2_CAMBER_BETA1', 'S3_CAMBER_BETA1',
             'S1_CAMBER_BETA2', 'S2_CAMBER_BETA2', 'S3_CAMBER_BETA2'],
            ['rb'] * 6
        ]
        index = pd.MultiIndex.from_tuples(
            list(zip(*arrays)), names=self.names
        )
        angles_throat = angles_throat.append(pd.DataFrame(beta_ab, index=index))

        arrays = [
            ['S1_THROAT_WIDTH', 'S2_THROAT_WIDTH', 'S3_THROAT_WIDTH'] * 2,
            ['gv'] * 3 + ['rb'] * 3
        ]
        index = pd.MultiIndex.from_tuples(
            list(zip(*arrays)), names=self.names
        )
        angles_throat = angles_throat.append(pd.DataFrame(self.throat_width, index=index))

        return angles_throat

    def _gamma_compute(self, ksi_rb=0.3, ksi_gv=0.3):

        gamma = {}
        r_te = {}
        r_le = {}

        for stage in list(range(1, self.stages_num + 1)):
            gamma[stage] = (
                90 - (self.alfa[stage][3] + ksi_gv *
                      (180 - (self.alfa[stage][0] + self.alfa[stage][3]))),

                90 - (self.alfa[stage][4] + ksi_gv *
                      (180 - (self.alfa[stage][1] + self.alfa[stage][4]))),

                90 - (self.alfa[stage][5] + ksi_gv *
                      (180 - (self.alfa[stage][2] + self.alfa[stage][5]))),

                self.beta[stage][3] - 90 + ksi_rb *
                (180 - (self.beta[stage][0] + self.beta[stage][3])),

                self.beta[stage][4] - 90 + ksi_rb *
                (180 - (self.beta[stage][1] + self.beta[stage][4])),

                self.beta[stage][5] - 90 + ksi_rb *
                (180 - (self.beta[stage][2] + self.beta[stage][5]))
            )
            r_le[stage] = (
                0.06 * (self.z[stage][1] - self.z[stage][0]) / math.cos(math.radians(gamma[stage][0])) / 2,
                0.06 * (self.z[stage][1] - self.z[stage][0]) / math.cos(math.radians(gamma[stage][1])) / 2,
                0.06 * (self.z[stage][1] - self.z[stage][0]) / math.cos(math.radians(gamma[stage][2])) / 2,
                0.06 * (self.z[stage][3] - self.z[stage][2]) / math.cos(math.radians(gamma[stage][3])) / 2,
                0.06 * (self.z[stage][3] - self.z[stage][2]) / math.cos(math.radians(gamma[stage][4])) / 2,
                0.06 * (self.z[stage][3] - self.z[stage][2]) / math.cos(math.radians(gamma[stage][5])) / 2,
            )
            r_te[stage] = (
                0.015 * (self.z[stage][1] - self.z[stage][0]) / math.cos(math.radians(gamma[stage][0])) / 2,
                0.015 * (self.z[stage][1] - self.z[stage][0]) / math.cos(math.radians(gamma[stage][1])) / 2,
                0.015 * (self.z[stage][1] - self.z[stage][0]) / math.cos(math.radians(gamma[stage][2])) / 2,
                0.015 * (self.z[stage][3] - self.z[stage][2]) / math.cos(math.radians(gamma[stage][3])) / 2,
                0.015 * (self.z[stage][3] - self.z[stage][2]) / math.cos(math.radians(gamma[stage][4])) / 2,
                0.015 * (self.z[stage][3] - self.z[stage][2]) / math.cos(math.radians(gamma[stage][5])) / 2,
            )

        arrays = [
            ['S1_CAMBER_GAMMA', 'S2_CAMBER_GAMMA', 'S3_CAMBER_GAMMA'] * 2,
            ['gv'] * 3 + ['rb'] * 3
        ]
        index = pd.MultiIndex.from_tuples(
            list(zip(*arrays)), names=self.names
        )
        gamma_df = pd.DataFrame(gamma, index=index)

        arrays = [
            ['S1_LE_RADIUS', 'S2_LE_RADIUS', 'S3_LE_RADIUS']*2,
            ['gv']*3 + ['rb']*3
        ]
        index = pd.MultiIndex.from_tuples(
            list(zip(*arrays)), names=self.names
        )
        gamma_df = gamma_df.append(pd.DataFrame(r_le, index=index))

        arrays = [
            ['S1_TE_RADIUS', 'S2_TE_RADIUS', 'S3_TE_RADIUS'] * 2,
            ['gv'] * 3 + ['rb'] * 3
        ]
        index = pd.MultiIndex.from_tuples(
            list(zip(*arrays)), names=self.names
        )
        gamma_df = gamma_df.append(pd.DataFrame(r_te, index=index))

        return gamma_df

    def number_of_blades(self):
        df = pd.DataFrame(self._data_formation())
        num_blades = {}

        for stage in list(range(1, self.stages_num + 1)):
            num_blades[stage] = df[stage]['CA.ZЛ']
        arrays = [
            ['NB'], ['gv']
        ]
        index = pd.MultiIndex.from_tuples(
            list(zip(*arrays)), names=self.names
        )
        z = pd.DataFrame(num_blades, index=index)

        for stage in list(range(1, self.stages_num + 1)):
            num_blades[stage] = df[stage]['PK.ZЛ']

        arrays = [
            ['NB'], ['rb']
        ]
        index = pd.MultiIndex.from_tuples(
            list(zip(*arrays)), names=self.names
        )
        z = z.append(
            pd.DataFrame(num_blades, index=index)
        )

        return z

    def generate_param_file(self, **kwargs):

        param_dir = os.path.join(self.root_dir, 'param')
        self.dir_handle.create_dir(param_dir)

        self.ab_data = self.ab_data.append(self._coordinates())
        self.ab_data = self.ab_data.append(self._twist_law())
        self.ab_data = self.ab_data.append(self._gamma_compute())
        self.ab_data = self.ab_data.append(self.number_of_blades())

        print(self.ab_data)

        for stage in self.ab_data.columns:
            for each in self.ab_data[stage].groupby('rs'):
                num_param = len(each[1]) + len(kwargs)
                param_str = ''
                param_str += f'# {num_param} parameters\n'
                param_str += '# Name\tExpression\tValue\n'
                # param_str += 'S1_SPAN\t0\nS2_SPAN\t0.5\nS3_SPAN\t1\t\n'
                file_name = os.path.join(param_dir, f'{each[0]}{stage}.param')
                for k, v in each[1].items():
                    param_str += f'{k[0]}\t{v}\t\n'
                if kwargs:
                    for k, v in kwargs.items():
                        param_str += f'{k}\t{v}\t\n'
                with open(file_name, 'w') as f:
                    f.write(param_str)

    def save_curve_files(self):
        curves_dir = os.path.join(self.root_dir, 'curves')
        log_msg, is_created = self.dir_handle.create_dir(curves_dir)
        if is_created:
            out_file_name_shroud = os.path.splitext(os.path.split(self.in_file)[1])[0] + '_shroud.dat'
            out_file_name_hub = os.path.splitext(os.path.split(self.in_file)[1])[0] + '_hub.dat'
            out_file_name_shroud = os.path.join(curves_dir, out_file_name_shroud)
            out_file_name_hub = os.path.join(curves_dir, out_file_name_hub)
