#!!!TODO проверить как это работает для разного количества ступеней, есть большое подозрение, что не будет

import pandas as pd


def Cycle(Data_for_cycle, df, range_, step):
    #Цикл, который разделяет строки на значени и параметры
    for i in range_:
        Param = df.loc[i].item().split()
        Value = df.loc[i + step].item().split()
        for j in range(len(Param)):
            Data_for_cycle[Param[j]] = float(Value[j])
    return Data_for_cycle


def data_for_stage(df, N):
    Data = {}
    # Блок расчета общих данных ступени
    low_boun = df.loc[df[0] == f'   {N}  CTУПEHЬ'].index[2] + 1
    high_boun = df.loc[df[0] == ' PEШETKA CA'].index[N-1] - 1
    Data = Cycle(Data, df, range_= range(low_boun,  high_boun, 2), step=1)

    # Блок расчета данных "Решетка СА", добавил еще радиусы скруглений и ширину сюда же
    CA = {}

    low_boun_CA_init = df.loc[df[0] == ' PEШETKA CA И PK'].index[N - 1] + 1
    high_boun_CA_init = low_boun_CA_init + 2
    Data['CA'] = Cycle(CA, df, range_=range(low_boun_CA_init, high_boun_CA_init, 2), step=1)

    low_boun_CA = df.loc[df[0] == ' PEШETKA CA'].index[N - 1] + 1
    high_boun_CA = df.loc[df[0] == ' PEШETKA PK'].index[N - 1] - 1
    Data['CA'] = Cycle(CA, df, range_=range(low_boun_CA, high_boun_CA, 2), step=1)


    # Блок расчета данных "Решетка РК", добавил еще радиусы скруглений и ширину сюда же
    PK = {}

    Data['PK'] = Cycle(PK, df, range_=range(low_boun_CA_init, high_boun_CA_init, 2), step=2)

    low_boun_PK = df.loc[df[0] == ' PEШETKA PK'].index[N - 1] + 1
    high_boun_PK = low_boun_PK + 10
    Data['PK'] = Cycle(PK, df, range_=range(low_boun_PK, high_boun_PK, 2), step=1)

    # Блок расчета начальных данных (тут только диаметры)
    Init = {}
    low_boun_init = df.loc[df[0] == f'   {N}  CTУПEHЬ'].index[0] + 1
    high_boun_init = low_boun_init + 2
    Data['Init'] = Cycle(Init, df, range_=range(low_boun_init, high_boun_init, 2), step=1)

    return Data


# df = pd.read_csv('St_base_m90.REZ', header=None, encoding='cp1251')

def number_of_stages(data_frame=None):
    l = data_frame.loc[data_frame[0] == '     Z    TAYMAX      Ш          N         DN       NK      ALFA2    ALFA2'].index[0]
    param = data_frame.loc[l + 1].item().split()
    return int(param[0])

DATA = {}

# for i in range(number_of_stages()):
#     DATA['Ступень {}'.format(i + 1)] = data_for_stage(sf, i + 1)


def value(stupen, parameter, *koleco):
    if koleco:
        return DATA[stupen][koleco[0]][parameter]
    else:
        return DATA[stupen][parameter]


for key in DATA.keys():
    print(f'{key}: {DATA[key]}')

# print(value('Ступень 3', 'D2', 'Init'))
# DATA['Параметры турбины'] = Cycle(DATA,df,range_=range(221,225,2))

# print(DATA['Ступень 3']['DE2'])
