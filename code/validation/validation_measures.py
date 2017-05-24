import random
import numpy as np
from sklearn.metrics import *
from collections import defaultdict

__author__ = 'Riccardo Guidotti'


def delta_k(real, pred):
    return len(set(pred)) - len(set(real))


def purity(real, pred):
    purity_val = 0
    cluster_count = defaultdict(lambda: [0] * len(set(real)))

    for p, r in zip(pred, real):
        cluster_count[p][r] += 1

    for label in cluster_count:
        purity_val += np.max(cluster_count[label])

    return 1.0 * purity_val / len(real)


def entropy(baskets):
    item_baskets = defaultdict(int)
    for basket in baskets:
        for item in basket:
            item_baskets[item] += 1
    entropy_val = 0
    size_cluster = len(baskets)
    for item in item_baskets:
        p = 1.0 * item_baskets[item] / size_cluster
        if p == 0 or p == 1.0:
            entropy_item = 0.0
        else:
            entropy_item = - ((p * np.log2(p)) + ((1-p) * np.log2(1-p)))
        entropy_val += entropy_item
    return entropy_val


def weighted_entropy(baskets_list, N):
    weighted_entropy_val = 0
    for baskets in baskets_list:
        entropy_val = entropy(baskets)
        n = len(baskets)
        weighted_entropy_val += 1.0 * n / N * entropy_val
    return weighted_entropy_val


def lisr(baskets, theta=0.5):
    item_baskets = defaultdict(int)
    for basket in baskets:
        for item in basket:
            item_baskets[item] += 1
    large_items = 0
    all_items = np.sum(item_baskets.values())
    for item in item_baskets:
        if item_baskets[item] >= len(baskets) * theta:
            large_items += item_baskets[item]
    return 1.0 * large_items / all_items


def weighted_lisr(baskets_list, N, theta=0.5):
    weighted_lisr_val = 0
    for baskets in baskets_list:
        lisr_val = lisr(baskets, theta)
        n = len(baskets)
        weighted_lisr_val += 1.0 * n / N * lisr_val
    return weighted_lisr_val


def wcd(baskets):
    item_baskets = defaultdict(int)
    for basket in baskets:
        for item in basket:
            item_baskets[item] += 1
    sum_occ = 0
    s_k = 0
    for item in item_baskets:
        sum_occ += item_baskets[item]**2
        s_k += item_baskets[item]
    n_k = len(baskets)
    if n_k == 0:
        return 0.0
    wcd_val = 1.0 * sum_occ / (s_k * n_k)
    return wcd_val


def weighted_wcd(baskets_list, N):
    weighted_wcd_val = 0
    for baskets in baskets_list:
        wcd_val = wcd(baskets)
        weighted_wcd_val += 1.0 * wcd_val * len(baskets) / N
    return weighted_wcd_val


def quality_practical(baskets, N, item_all_baskets):
    item_baskets = defaultdict(int)
    for basket in baskets:
        for item in basket:
            item_baskets[item] += 1
    quality_val = 0
    size_cluster = len(baskets)
    if size_cluster == 0:
        return 0.0
    for item in item_baskets:
        items_in_cluster = item_baskets[item]
        items_in_db = item_all_baskets[item]
        items_not_in_cluster = items_in_db - items_in_cluster
        w_item_cluster = 1.0 * (1.0 * items_in_cluster / size_cluster) * (1.0 * items_in_cluster / (items_in_cluster + items_not_in_cluster))
        w_item_db = 1.0 * items_in_db * (N - items_in_db + 1) / N
        quality_val += (1.0 * items_in_cluster * w_item_cluster * w_item_db)
    quality_val /= size_cluster

    return quality_val


def weighted_quality_practical(baskets_list, N):
    item_all_baskets = defaultdict(int)
    for baskets in baskets_list:
        for basket in baskets:
            for item in basket:
                item_all_baskets[item] += 1
    weighted_quality_val = 0
    for baskets in baskets_list:
        quality_val = quality_practical(baskets, N, item_all_baskets)
        n = len(baskets)
        weighted_quality_val += 1.0 * n / N * quality_val
    return weighted_quality_val


def quality_atdc(baskets, N, item_all_baskets):
    item_baskets = defaultdict(int)
    for basket in baskets:
        for item in basket:
            item_baskets[item] += 1
    quality_val = 0
    n = len(baskets)
    if n == 0:
        return 0.0
    for item in item_baskets:
        n_a = item_baskets[item]
        N_a = item_all_baskets[item]
        quality_val += ((1.0 * n_a / n)** 2 - (1.0 * N_a / N)**2)

    return quality_val


def weighted_quality_atdc(baskets_list, N):
    item_all_baskets = defaultdict(int)
    for baskets in baskets_list:
        for basket in baskets:
            for item in basket:
                item_all_baskets[item] += 1
    weighted_quality_val = 0
    for baskets in baskets_list:
        quality_val = quality_atdc(baskets, N, item_all_baskets)
        n = len(baskets)
        weighted_quality_val += 1.0 * n / N * quality_val
    return weighted_quality_val


def jaccard_bitarray(a, b):
    intersection_cardinality = (a & b).count()
    union_cardinality = (a | b).count()
    if union_cardinality == 0:
        return 1.0
    else:
        return 1.0 - 1.0 * intersection_cardinality / union_cardinality


def silhouette(baskets, labels, N, nsample=float('infinity'), niter=10):

    if nsample > N:
        baskets_keys_sample = baskets.keys()
        labels_sample = np.asarray(labels)
        distances = list()
        for b1 in baskets_keys_sample:
            distance = list()
            for b2 in baskets_keys_sample:
                dist = jaccard_bitarray(baskets[b1], baskets[b2])
                distance.append(dist)
            distances.append(np.asarray(distance))

        distances = np.asarray(distances)

        if 1 < len(set(labels_sample)) < N:
            return silhouette_score(distances, labels_sample, metric='precomputed')
        else:
            return 0.0

    s_list = list()
    for i in range(0, niter):
        baskets_keys_sample = random.sample(baskets.keys(), nsample)
        labels_sample = np.asarray(map(labels.__getitem__, baskets_keys_sample))

        distances = list()
        for b1 in baskets_keys_sample:
            distance = list()
            for b2 in baskets_keys_sample:
                dist = jaccard_bitarray(baskets[b1], baskets[b2])
                distance.append(dist)
            distances.append(np.asarray(distance))

        distances = np.asarray(distances)

        if 1 < len(set(labels_sample)) < nsample:
            s_val = silhouette_score(distances, labels_sample, metric='precomputed')
        else:
            s_val = 0.0
        s_list.append(s_val)

    return np.mean(s_list)

