from scipy import stats
from bitarray import bitarray
from collections import defaultdict

__author__ = 'Riccardo Guidotti'


def jaccard_bitarray(a, b):
    intersection_cardinality = (a & b).count()
    union_cardinality = (a | b).count()
    if union_cardinality == 0:
        return 1.0
    else:
        return 1.0 - 1.0 * intersection_cardinality / union_cardinality


def sample_size(population, conf_interval, conf_level=0.95, prob=0.5):
    if conf_interval < 0.01:
        return population

    zscore = stats.norm.ppf(1-(1-conf_level)/2)
    ss = ((zscore ** 2) * prob * (1-prob)) / (conf_interval ** 2)
    new_ss = ss / (1 + (ss - 1) / population)

    return int(round(new_ss))


def calculate_item_baskets(baskets, nbaskets):
    item_baskets = defaultdict(lambda: nbaskets * bitarray('0'))
    for b in baskets:
        for item in range(0, len(baskets[b])):
            if baskets[b][item]:
                item_baskets[item][b] = 1
    return item_baskets


def count_items(baskets):

    item_set = set()
    for b in baskets:
        for item in range(0, len(baskets[b])):
            if baskets[b][item]:
                item_set.add(item)

    return len(item_set)


def remap_items(baskets_list):
    map_item_newitem = dict()
    map_newitem_item = dict()
    baskets_list_new = list()
    for basket_old in baskets_list:
        basket_new = list()
        for item in basket_old:
            if item not in map_item_newitem:
                newitem = len(map_item_newitem)
                map_item_newitem[item] = newitem
                map_newitem_item[newitem] = item
            newitem = map_item_newitem[item]
            basket_new.append(newitem)

        baskets_list_new.append(basket_new)
    return baskets_list_new, map_newitem_item, map_item_newitem


def remap_items_back(baskets_list, map_newitem_item):

    baskets_list_old = list()
    for basket_new in baskets_list:
        basket_old = list()
        for newitem in basket_new:
            item = map_newitem_item[newitem]
            basket_old.append(item)

        baskets_list_old.append(basket_old)
    return baskets_list_old


def basket_bitarray_to_list(baskets_bitarray):

    baskets_list = dict()

    for b in baskets_bitarray:
        baskets_list[b] = list()
        for item in range(0, len(baskets_bitarray[b])):
            if baskets_bitarray[b][item]:
                baskets_list[b].append(item)

    return baskets_list


def basket_list_to_bitarray(baskets_list, nitems):

    baskets_bitarray = dict()
    for b, basket in enumerate(baskets_list):
        baskets_bitarray[b] = nitems * bitarray('0')
        for item in baskets_list[b]:
            baskets_bitarray[b][item] = 1

    return baskets_bitarray
