# -*- coding: utf-8 -*-
import pip
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

    def __init__(self, in_file=None, root_dir=os.getcwd(), parameters=None, rs_name=None, coef=1, **laws):

        self.dir_handle = CommonClass()
        self.root_dir = root_dir
        self.names = ['param', 'rs']
        self.in_file = in_file
        self.parameters = parameters
        self.rs_name = rs_name
        self.coef = coef
        self.hub = {}
        self.tip = {}
        self.b1 = {}
        self.z = {}
        self.ref_trace = {}
        self.r_tip = {}
        self.r_hub = {}
        self.beta = {}
        self.alfa = {}
        self.throat_width = {}
        self.stages_num = 0
        self.data_frame = pd.read_csv(in_file, header=None, encoding="cp1251")
        self.data = {}
        self.ab_data = pd.DataFrame()
        self.laws = laws

    @staticmethod
    def line_equation(x1, y1, x2, y2, yi):
        xi = x2 - ((y2 - yi) / (y2 - y1)) * (x2 - x1)
        return xi

    @staticmethod
    def slope_angle(x1, y1, x2, y2):
        slope = math.degrees(math.atan((y2 - y1) / (x2 - x1)))
        return slope

    @staticmethod
    def get_conicity(d_rel, temperature):
        if temperature > 473.15:
            a, b = 0.0847, 0.2202
        else:
            a, b = 0.0918, -0.01#-0.0513
        conicity = a * d_rel + b
        if conicity > 1:
            conicity = 1
        return conicity

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

    def _coordinates(self, dz_=1.2, sweep_beta=0, stacking_law=0, surface_setup=0):
        if 'stacking_law' in self.laws.keys():
            stacking_law = self.laws['stacking_law']
        if 'surface_setup' in self.laws.keys():
            surface_setup = self.laws['surface_setup']

        # Stacking laws:
        # 0 - Center of gravity
        # 1 - Leading edge
        # 2 - Trailing edge
        # in any other case - Trailing edge

        stacking_point = {
            0: 'GRAVITY_CENTER',
            1: 'LEADING_EDGE',
            2: 'TRAILING_EDGE'
        }
        stack = {}

        df = pd.DataFrame(self._data_formation())

        d_throat_tip = {}
        d_throat_hub = {}
        z0 = 0
        dz = {}
        z_coordinates = {}
        b_ax = {}

        zth = {}
        r_tip = {}
        r_hub = {}

        d_tip = {}
        d_hub = {}
        sweep = {}

        for stage in list(range(1, self.stages_num + 1)):

            conicity = self.get_conicity(
                df[stage]['Init.D2'] / df[stage]['PK.HЛ'],
                df[stage]['T2*']
            )
            b_mid = df[stage]['PK.L'] * self.coef * math.cos(math.radians(df[stage]['CA.ГAMMA']))
            b_ax[stage] = (
                2 * b_mid / (1 + conicity),
                b_mid,
                conicity * 2 * b_mid / (1 + conicity)
            )
            d_tip[stage] = []
            d_hub[stage] = []
            d_throat_tip[stage] = []
            d_throat_hub[stage] = []

            self.z[stage] = []

            dz[stage] = [
                dz_ * df[stage]['CA.T'] * self.coef *
                math.tan(math.radians(df[stage]['AL1'])),
                dz_ * df[stage]['PK.T'] * self.coef *
                math.tan(math.radians(df[stage]['BE2']))
            ]

            if stage > 1:
                z0 += (df[stage]['CA.L'] * self.coef * math.cos(math.radians(df[stage]['CA.ГAMMA']))
                       + dz[stage][0] + df[stage]['PK.L'] * self.coef *
                       math.cos(math.radians(df[stage]['PK.ГAMMA'])) + dz[stage][1])

            throat_gv = (0.5 * df[stage]['CA.A/T'] * df[stage]['CA.T'] * self.coef *
                         math.cos(math.radians(df[stage]['AL1'])))
            throat_rb = 0.5 * df[stage]['PK.A/T'] * df[stage]['PK.T'] * \
                        math.cos(math.radians(df[stage]['BE2']))

            z_throat_gv = z0 + df[stage]['CA.L'] * self.coef * \
                          math.sin(math.radians(df[stage]['CA.ГAMMA'])) - throat_gv
            z1 = z0 + df[stage]['CA.L'] * self.coef * math.sin(math.radians(df[stage]['CA.ГAMMA']))
            z2 = z1 + dz[stage][0]
            z_throat_rb = z2 + df[stage]['PK.L'] * self.coef * \
                          math.sin(math.radians(df[stage]['PK.ГAMMA'])) - throat_rb
            z3 = z2 + df[stage]['PK.L'] * self.coef * math.sin(math.radians(df[stage]['PK.ГAMMA']))

            z_coordinates[stage] = [z0, z_throat_gv, z1, z2, z_throat_rb, z3]

            d_throat_tip[stage].append((df[stage]['Init.D1'] + df[stage]['CA.HЛ']) * self.coef)
            d_throat_tip[stage].append((df[stage]['Init.D2'] + df[stage]['PK.HЛ']) * self.coef)
            d_throat_hub[stage].append((df[stage]['Init.D1'] - df[stage]['CA.HЛ']) * self.coef)
            d_throat_hub[stage].append((df[stage]['Init.D2'] - df[stage]['PK.HЛ']) * self.coef)

            self.z[stage] = tuple([z0, z1, z2, z3])

            zth[stage] = tuple([
                -0.005 if stage == 1 else z0 - z0 * 0.1,
                z1 + z1 * 0.1, z2 - z2 * 0.1, z3 + z3 * 0.1])

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
            if surface_setup == 0:
                self.ref_trace[stage] = tuple([
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
            elif surface_setup == 1:
                self.ref_trace[stage] = tuple([
                    self.z[stage][0], self.z[stage][0],
                    self.r_tip[stage][0] if self.r_tip[stage][0] > self.r_tip[stage][1]
                    else self.r_tip[stage][1],
                    self.r_hub[stage][0] if self.r_hub[stage][0] < self.r_hub[stage][1]
                    else self.r_hub[stage][1],
                    0, 0,
                    self.z[stage][2], self.z[stage][2],
                    self.r_tip[stage][2] if self.r_tip[stage][2] > self.r_tip[stage][3]
                    else self.r_tip[stage][3],
                    self.r_hub[stage][2] if self.r_hub[stage][2] < self.r_hub[stage][3]
                    else self.r_hub[stage][3],
                    0, 0
                ])
            else:
                pass

            stack[stage] = [
                stacking_point[stacking_law[0]],
                stacking_point[stacking_law[1]]
            ]
            if stacking_law[0] == 0:
                sweep[stage] = [self.z[stage][0] + (self.z[stage][1] - self.z[stage][0]) / 2,
                                sweep_beta]
            elif stacking_law[0] == 1:
                sweep[stage] = [self.z[stage][0], sweep_beta]
            elif stacking_law[0] == 2:
                sweep[stage] = [self.z[stage][1], sweep_beta]

            if stacking_law[1] == 0:
                sweep[stage].append(self.z[stage][2] + (self.z[stage][3] - self.z[stage][2]) / 2)
                sweep[stage].append(sweep_beta)
            elif stacking_law[1] == 1:
                sweep[stage].append(self.z[stage][2])
                sweep[stage].append(sweep_beta)
            elif stacking_law[1] == 2:
                sweep[stage].append(self.z[stage][1])
                sweep[stage].append(sweep_beta)

            self.b1[stage] = (
                self.z[stage][1] - self.z[stage][0],
                self.z[stage][1] - self.z[stage][0],
                self.z[stage][1] - self.z[stage][0],
                b_ax[stage][0], b_ax[stage][1], b_ax[stage][2]
                # self.z[stage][3] - self.z[stage][2],
                # self.z[stage][3] - self.z[stage][2],
                # self.z[stage][3] - self.z[stage][2]
            )

        arrays = [
            ['HUB_R1', 'HUB_R2'] * 2,
            ['gv', 'gv', 'rb', 'rb']
        ]
        index = pd.MultiIndex.from_tuples(
            list(zip(*arrays)), names=self.names
        )
        coordinates = pd.DataFrame(r_hub, index=index)

        arrays = [
            ['SHROUD_R1', 'SHROUD_R2'] * 2,
            ['gv', 'gv', 'rb', 'rb']
        ]
        index = pd.MultiIndex.from_tuples(
            list(zip(*arrays)), names=self.names
        )
        coordinates = coordinates.append(pd.DataFrame(r_tip, index=index))

        arrays = [
            ['HUB_Z1', 'HUB_Z2'] * 2,
            ['gv', 'gv', 'rb', 'rb']
        ]
        index = pd.MultiIndex.from_tuples(
            list(zip(*arrays)), names=self.names
        )
        coordinates = coordinates.append(pd.DataFrame(zth, index=index))

        arrays = [
            ['SHROUD_Z1', 'SHROUD_Z2'] * 2,
            ['gv', 'gv', 'rb', 'rb']
        ]
        index = pd.MultiIndex.from_tuples(
            list(zip(*arrays)), names=self.names
        )
        coordinates = coordinates.append(pd.DataFrame(zth, index=index))

        # if surface_setup == 0:
        arrays = [
            ['REF_TRACE_TIP_Z', 'REF_TRACE_HUB_Z',
             'REF_TRACE_TIP_R', 'REF_TRACE_HUB_R',
             'REF_TRACE_TIP_ALPHA', 'REF_TRACE_HUB_ALPHA'] * 2,
            ['gv'] * 6 + ['rb'] * 6
        ]
        index = pd.MultiIndex.from_tuples(
            list(zip(*arrays)), names=self.names
        )
        coordinates = coordinates.append(pd.DataFrame(self.ref_trace, index=index))

        arrays = [
            ['Z_STACKING_SWEEP', 'SWEEP_BETA'] * 2,
            ['gv'] * 2 + ['rb'] * 2
        ]
        index = pd.MultiIndex.from_tuples(
            list(zip(*arrays)), names=self.names
        )
        coordinates = coordinates.append(pd.DataFrame(sweep, index=index))

        arrays = [
            ['STACKING_POINT'] * 2,
            ['gv'] + ['rb']
        ]
        index = pd.MultiIndex.from_tuples(
            list(zip(*arrays)), names=self.names
        )
        coordinates = coordinates.append(pd.DataFrame(stack, index=index))
        arrays = [
            ['S1_DZ_PRIM', 'S2_DZ_PRIM', 'S3_DZ_PRIM'] * 2,
            ['gv'] * 3 + ['rb'] * 3
        ]

        index = pd.MultiIndex.from_tuples(
            list(zip(*arrays)), names=self.names
        )

        coordinates = coordinates.append(pd.DataFrame(self.b1, index=index))

        self.save_curve_files(z_coordinates, d_tip, d_hub)
        return coordinates

    def _twist_law(self, twist_law=0, profile_law=0):

        if 'profile_law' in self.laws.keys():
            profile_law = self.laws['profile_law']
        if 'twist_law' in self.laws.keys():
            twist_law = self.laws['twist_law']

        # twist_law: 0 - twisted by constant circulation law;
        # twist_law: 1 - twisted by constant reactivity
        # profile_law: 0 - Throat base law
        # profile_law: 1 - Legacy law

        df = pd.DataFrame(self._data_formation())

        alfa_2 = {}
        for stage in list(range(1, self.stages_num + 1)):

            alfa2_mid = float(df[stage]['AL2'])
            alfa1_mid = float(df[stage]['AL1'])
            alfa0_mid = float(df[stage]['AL0'])
            h0_tot = df[stage]['HAД'] * 9.80665  # J/kg
            hu = df[stage]['LU'] * 9.80665  # J/kg
            ro_mid = df[stage]['PO']
            ro_hub = df[stage]['POK']
            phi = df[stage]['FИ']
            u1_mid = df[stage]['V1']

            # [s0, s1, s2, s3]
            # [hub, mid, shroud]
            r_mid = [
                ((self.r_tip[stage][0] + self.r_hub[stage][0]) / 2) / self.coef,
                ((self.r_tip[stage][1] + self.r_hub[stage][1]) / 2) / self.coef,
                ((self.r_tip[stage][2] + self.r_hub[stage][2]) / 2) / self.coef,
                ((self.r_tip[stage][3] + self.r_hub[stage][3]) / 2) / self.coef
            ]
            omega = u1_mid / r_mid[2]

            if twist_law == 0:
                alfa_1 = [
                    math.degrees(math.atan(((self.r_hub[stage][1] / self.coef) / r_mid[1]) *
                                           math.tan(math.radians(alfa1_mid)))),
                    alfa1_mid,
                    math.degrees(math.atan(((self.r_tip[stage][1] / self.coef) / r_mid[1]) *
                                           math.tan(math.radians(alfa1_mid))))
                ]
                ro = [
                    ro_hub,
                    ro_mid,
                    1 - (1 - ro_mid) * (r_mid[1] / (self.r_tip[stage][1] / self.coef)) ** 2,
                ]
                h01 = [h0_tot * (1 - ro_i) for ro_i in ro]
                c1 = [phi * (2 * h01_i) ** 0.5 for h01_i in h01]
                c1z = c1[1] * math.sin(math.radians(alfa1_mid))
                c1u = [a[0] * math.cos(math.radians(a[1])) for a in zip(c1, alfa_1)]
                u1 = [
                    u1_mid * ((self.r_hub[stage][2] / self.coef) / r_mid[2]),
                    u1_mid,
                    u1_mid * ((self.r_tip[stage][2] / self.coef) / r_mid[2])
                ]
                u2 = [
                    omega * self.r_hub[stage][3] / self.coef,
                    omega * r_mid[3],
                    omega * self.r_tip[stage][3] / self.coef
                ]

                beta_1 = [math.degrees(math.atan(c1z / (a[0] - a[1]))) for a in zip(c1u, u1)]
                beta_1 = [beta_1_i if beta_1_i > 0 else beta_1_i + 180 for beta_1_i in beta_1]
                c2u = [(hu - (a[0] * a[1])) / a[2] for a in zip(u1, c1u, u2)]
                c2z = c2u[1] * math.tan(math.radians(alfa2_mid))
                beta_2 = [math.degrees(math.atan(c2z / (a[0] + a[1]))) for a in zip(c2u, u2)]
                alfa_2[stage] = [math.degrees(math.atan(c2z / c2u_i)) for c2u_i in c2u]

                if stage == 1:
                    alfa_0 = [alfa0_mid] * 3
                else:
                    alfa_0 = alfa_2[stage - 1]

                self.beta[stage] = beta_1 + beta_2
                self.alfa[stage] = alfa_0 + alfa_1

            elif twist_law == 1:
                n = (math.cos(math.radians(alfa1_mid))) ** 2
                # [hub, mid, shroud]
                ro = [
                    ro_hub,
                    ro_mid,
                    1 - (1 - ro_mid) * (r_mid[1] / (self.r_tip[stage][1] / self.coef)) ** n,
                ]
                h01 = [h0_tot * (1 - ro_i) for ro_i in ro]
                c1 = [phi * (2 * h01_i) ** 0.5 for h01_i in h01]
                c1z = [c1_i * math.sin(math.radians(alfa1_mid)) for c1_i in c1]
                c1u = [c1_i * math.cos(math.radians(alfa1_mid)) for c1_i in c1]

                u1 = [
                    u1_mid * (self.r_hub[stage][2] / self.coef) / r_mid[2],
                    u1_mid,
                    u1_mid * (self.r_tip[stage][2] / self.coef) / r_mid[2]
                ]
                u2 = [
                    omega * self.r_hub[stage][3] / self.coef,
                    omega * r_mid[3],
                    omega * self.r_tip[stage][3] / self.coef
                ]
                beta_1 = [math.degrees(math.atan(a[0] / (a[1] - a[2]))) for a in zip(c1z, c1u, u1)]
                beta_1 = [beta_1_i if beta_1_i > 0 else beta_1_i + 180 for beta_1_i in beta_1]
                alfa_1 = [alfa1_mid] * 3

                c2u = [(hu - a[0] * a[1]) / a[2] for a in zip(u1, c1u, u2)]
                c2_mid = c2u[1] / math.cos(math.radians(alfa2_mid))
                n1 = (n - 1) / n
                n2 = 2 * n
                n3 = 1 / n
                n4 = (n - 1) / (n + 1)
                n5 = 2 / (n + 1)
                n6 = n + 1
                r1 = r_mid[3] / (self.r_hub[stage][3] / self.coef)
                r2 = r_mid[3] / (self.r_tip[stage][3] / self.coef)
                dc = (c1u[1] - c2u[1])
                c2z = [
                    (c2_mid ** 2 - (c1u[1] ** 2) * (n1 * r1 ** n2 + n3) +
                     2 * c1u[1] * dc * (n4 * r1 ** n6 + n5) - dc ** 2) ** 0.5,
                    (c2_mid ** 2 - c2u[1] ** 2) ** 0.5,
                    (c2_mid ** 2 - (c1u[1] ** 2) * (n1 * r2 ** n2 + n3) +
                     2 * c1u[1] * dc * (n4 * r2 ** n6 + n5) - dc ** 2) ** 0.5
                ]
                beta_2 = [math.degrees(math.atan(a[0] / (math.fabs(a[1]) + a[2]))) for a in zip(c2z, c2u, u2)]
                alfa_2[stage] = [math.degrees(math.atan(a[0] / math.fabs(a[1]))) for a in zip(c2z, c2u)]

                if stage == 1:
                    alfa_0 = [alfa0_mid] * 3
                else:
                    alfa_0 = alfa_2[stage - 1]

                self.beta[stage] = beta_1 + beta_2
                self.alfa[stage] = alfa_0 + alfa_1

            if profile_law == 0:
                # Throat width
                # throat_width[stage][GV(hub, mid, tip), RB(hub, mid, tip)]
                self.throat_width[stage] = [
                    math.sin(math.radians(alfa_1[0])) * 2 *
                    math.pi * self.r_hub[stage][1] / df[stage]['CA.ZЛ'],

                    math.sin(math.radians(alfa_1[1])) * 2 *
                    math.pi * r_mid[1] / df[stage]['CA.ZЛ'],

                    math.sin(math.radians(alfa_1[2])) * 2 *
                    math.pi * self.r_tip[stage][1] / df[stage]['CA.ZЛ'],

                    math.sin(math.radians(beta_2[0])) * 2 *
                    math.pi * self.r_hub[stage][3] / df[stage]['PK.ZЛ'],

                    math.sin(math.radians(beta_2[1])) * 2 *
                    math.pi * r_mid[3] / df[stage]['PK.ZЛ'],

                    math.sin(math.radians(beta_2[2])) * 2 *
                    math.pi * self.r_tip[stage][3] / df[stage]['PK.ZЛ']
                ]
            elif profile_law == 1:
                pass
            else:
                pass

        alfa_ab = {}
        for stage in self.alfa.keys():
            alfa_ab[stage] = [(90 - x) for x in self.alfa[stage][:3]]
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

        if profile_law == 0:
            arrays = [
                ['S1_THROAT_WIDTH', 'S2_THROAT_WIDTH', 'S3_THROAT_WIDTH'] * 2,
                ['gv'] * 3 + ['rb'] * 3
            ]
            index = pd.MultiIndex.from_tuples(
                list(zip(*arrays)), names=self.names
            )
            angles_throat = angles_throat.append(pd.DataFrame(self.throat_width, index=index))
        elif profile_law == 1:
            pass
        else:
            pass

        return angles_throat

    def _gamma_compute(self, ksi_rb=0.3, ksi_gv=0.3):

        df = pd.DataFrame(self._data_formation())

        gamma = {}
        r_te = {}
        r_le = {}
        for stage in list(range(1, self.stages_num + 1)):
            t_rb = [
                2 * self.r_hub[stage][2] * math.pi / df[stage]["PK.ZЛ"],
                df[stage]["PK.T/L"] * df[stage]["PK.L"] * self.coef,
                2 * self.r_tip[stage][2] * math.pi / df[stage]["PK.ZЛ"]
            ]

            b_rb = [a[0] / a[1] for a in zip(t_rb, [df[stage]["PK.T/L"]]*3)]

            gamma_rb = [a[0] / a[1] if a[0] / a[1] < 1 else 1
                        for a in zip(self.b1[stage][3:], b_rb)]


            gamma[stage] = (
                90 - (self.alfa[stage][3] + ksi_gv *
                      (180 - (self.alfa[stage][0] + self.alfa[stage][3]))),

                90 - (self.alfa[stage][4] + ksi_gv *
                      (180 - (self.alfa[stage][1] + self.alfa[stage][4]))),

                90 - (self.alfa[stage][5] + ksi_gv *
                      (180 - (self.alfa[stage][2] + self.alfa[stage][5]))),

                math.degrees(math.asin(gamma_rb[0])) - 90,
                math.degrees(math.asin(gamma_rb[1])) - 90,
                math.degrees(math.asin(gamma_rb[2])) - 90

                # self.beta[stage][3] - 90 + ksi_rb *
                # (180 - (self.beta[stage][0] + self.beta[stage][3])),
                #
                # self.beta[stage][4] - 90 + ksi_rb *
                # (180 - (self.beta[stage][1] + self.beta[stage][4])),
                #
                # self.beta[stage][5] - 90 + ksi_rb *
                # (180 - (self.beta[stage][2] + self.beta[stage][5]))
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
            ['S1_LE_RADIUS', 'S2_LE_RADIUS', 'S3_LE_RADIUS'] * 2,
            ['gv'] * 3 + ['rb'] * 3
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

    def generate_par_file(self, **kwargs):
        root = os.getcwd()
        par_dir = os.path.join(self.root_dir, 'par')
        par_sample = os.path.join(root, 'par_samples', 'sample.par')
        log_msg, is_created = self.dir_handle.create_dir(par_dir)

        self.ab_data = self.ab_data.append(self._coordinates())
        self.ab_data = self.ab_data.append(self._twist_law())
        self.ab_data = self.ab_data.append(self._gamma_compute())
        self.ab_data = self.ab_data.append(self.number_of_blades())

        if is_created:
            for stage in self.ab_data.columns:
                for each in self.ab_data[stage].groupby('rs'):
                    par_file = os.path.join(par_dir, f'{each[0]}{stage}.par')
                    try:
                        with open(par_sample, 'r') as sample:
                            is_param_block = False
                            param = False
                            with open(par_file, 'w') as par:
                                for line in sample:
                                    words = line.split()
                                    if len(words) == 2:
                                        if (str(words[0]) == "NI_BEGIN" and
                                                str(words[1]) in ["PARAMETERS", "PARAMETER", "STACKING_POINT"]):
                                            is_param_block = True
                                        elif (str(words[0]) == "NI_END" and
                                                str(words[1]) in ["PARAMETERS", "PARAMETER", "STACKING_POINT"]):
                                            is_param_block = False
                                            param = False
                                        elif is_param_block:
                                            if str(words[0]) in ["NAME"]:
                                                for k in each[1].keys():
                                                    if k[0] in words:
                                                        key = k[0]
                                                        param = True
                                            elif str(words[0]) in ["STACKING_POINT"]:
                                                for k in each[1].keys():
                                                    if k[0] in words:
                                                        key = k[0]
                                                        param = True
                                                line = line.replace(str(words[1]), str(each[1][key][0]))
                                            elif str(words[0]) in ["VALUE"]:
                                                if param:
                                                    line = line.replace(str(words[1]), str(each[1][key][0]))
                                            elif str(words[0]) == "LIMIT_MIN":
                                                if param:
                                                    line = line.replace(
                                                        str(words[1]), str(
                                                            each[1][key][0] - each[1][key][0] * 0.1
                                                        )
                                                    )
                                            elif str(words[0]) == "LIMIT_MAX":
                                                if param:
                                                    line = line.replace(
                                                        str(words[1]), str(
                                                            each[1][key][0] + each[1][key][0] * 0.1
                                                        )
                                                    )
                                            elif str(words[0]) == "VALUE_MIN":
                                                if param:
                                                    line = line.replace(
                                                        str(words[1]), str(
                                                            each[1][key][0] - each[1][key][0] * 0.1
                                                        )
                                                    )
                                            elif str(words[0]) == "VALUE_MAX":
                                                if param:
                                                    line = line.replace(
                                                        str(words[1]), str(
                                                            each[1][key][0] + each[1][key][0] * 0.1
                                                        )
                                                    )
                                    par.write(line)

                    except Exception as ex:
                        print(ex)

    def save_curve_files(self, *args):
        curves_dir = os.path.join(self.root_dir, 'curves')
        log_msg, is_created = self.dir_handle.create_dir(curves_dir)
        z = args[0]
        r_shroud = {key: [float(arg_i) / 2 for arg_i in arg] for key, arg in args[1].items()}
        r_hub = {key: [float(arg_i) / 2 for arg_i in arg] for key, arg in args[2].items()}

        if is_created:
            out_file_shroud = os.path.join(curves_dir, 'shroud.dat')
            out_file_hub = os.path.join(curves_dir, 'hub.dat')
            with open(out_file_shroud, 'w') as fsh:
                for stage in z.keys():
                    c = zip(z[stage], r_shroud[stage])
                    for e in c:
                        fsh.write(f'{e[0]},{e[1]},0\n')
            with open(out_file_hub, 'w') as fhub:
                for stage in z.keys():
                    c = zip(z[stage], r_hub[stage])
                    for e in c:
                        fhub.write(f'{e[0]},{e[1]},0\n')
