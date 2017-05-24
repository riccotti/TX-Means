import pandas as pd
import numpy as np
from scipy.stats import mode

__author__ = 'Riccardo Guidotti'

def read_synthetic_data(filename):

    baskets = list()

    data = open(filename, 'r')
    data.readline()
    for row in data:
        fields = row.rstrip().split(';')
        if len(fields[2]) > 0:
            baskets.append((map(int, fields[2].split(' ')), int(fields[1])))

    data.close()
    return baskets


def write_syntetic_data(filename, baskets):
    out = open(filename, 'w')
    out.write('ID;CLASS;EVENTS\n')
    for index in range(0, len(baskets)):
        basket = baskets[index]
        out.write('%d;%d;%s\n' % (index, basket[1], ' '.join(str(x) for x in basket[0])))

    out.flush()
    out.close()


def read_uci_data(filename, class_index=0, delimiter=',', missing_symbol='?', header=True, skipcolumnsindex=set()):
    df = pd.read_csv(filename, skipinitialspace=True)
    index_mode = dict()

    for k, index in zip(df.columns, range(0, len(df.columns))):
        df[k] = df[k].replace('?', np.nan)
        mode_value = mode(df[k])[0][0]
        df[k] = df[k].fillna(mode_value)
        index_mode[index] = mode_value

    baskets = list()

    data = open(filename, 'r')

    map_item_newitem = dict()
    map_newitem_item = dict()
    map_class_newclass = dict()
    map_newclass_class = dict()

    if header:
        data.readline()
    for row in data:
        categories = row.rstrip().split(delimiter)
        basket = list()
        basket_class = None
        for index in range(0, len(categories)):

            if index in skipcolumnsindex:
                continue

            if index == class_index:
                cclass = categories[index]
                if categories[index] not in map_class_newclass:
                    newclass = len(map_class_newclass)
                    map_class_newclass[cclass] = newclass
                    map_newclass_class[newclass] = cclass
                basket_class = map_class_newclass[cclass]
                continue

            if categories[index] == missing_symbol:
                categories[index] = index_mode[index]

            item = (index, categories[index])
            if item not in map_item_newitem:
                newitem = len(map_item_newitem)
                map_item_newitem[item] = newitem
                map_newitem_item[newitem] = item
            newitem = map_item_newitem[item]
            basket.append(newitem)

        if len(basket) > 0:
            baskets.append((basket, basket_class))

    data.close()

    maps = {
        'map_item_newitem': map_item_newitem,
        'map_newitem_item': map_newitem_item,
        'map_class_newclass': map_class_newclass,
        'map_newclass_class': map_newitem_item,
    }
    return baskets, maps
