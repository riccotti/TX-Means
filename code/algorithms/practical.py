import random
import datetime
from util import *

__author__ = 'Riccardo Guidotti'


def calculate_w_m_D(item_baskets, nbaskets, nitems):
    w_m_D = defaultdict(lambda: nitems * bitarray('0'))
    n = nbaskets
    for item in item_baskets:
        ntr_m_D = item_baskets[item].count()
        w_m_D[item] = 1.0 * (ntr_m_D * (n - ntr_m_D + 1)) / n
    return w_m_D


def w_m_C(cluster_info, item):
    fact1 = 1.0 * cluster_info['ntr_m_Cj'][item] / cluster_info['nj']
    fact2 = 1.0 * cluster_info['ntr_m_Cj'][item] / (cluster_info['ntr_m_Cj'][item] + cluster_info['ntr_m_D_Cj'][item])
    return fact1 * fact2


def calculate_quality_cluster(cluster_info, w_m_D, nbaskets, nitems):
    quality_sum = 0
    for item in range(0, nitems):
        val = (cluster_info['ntr_m_Cj'][item] + 1) * w_m_C(cluster_info, item) * w_m_D[item]
        quality_sum += val

    return quality_sum


def calculate_quality(quality_dict, nbaskets):
    return sum(quality_dict.values()) / nbaskets


def calculate_quality_dict(clusters_info, w_m_D, nbaskets, nitems):
    quality_dict = dict()
    for cid, cluster_info in clusters_info.iteritems():
        val = calculate_quality_cluster(cluster_info, w_m_D, nbaskets, nitems)
        quality_dict[cid] = val
    return quality_dict


def calculate_cluster_info_update(cluster_info, w_m_D, nbaskets, nitems, basket, update_factor=1):
    update_cluster_info = dict()
    update_cluster_info['nj'] = cluster_info['nj'] + update_factor
    if update_cluster_info['nj'] < 1:
        return None
    update_cluster_info['ntr_m_Cj'] = [0] * nitems
    update_cluster_info['ntr_m_D_Cj'] = [0] * nitems

    for item in range(0, nitems):
        if basket[item] == 1:
            update_cluster_info['ntr_m_Cj'][item] = cluster_info['ntr_m_Cj'][item] + update_factor
            update_cluster_info['ntr_m_D_Cj'][item] = cluster_info['ntr_m_D_Cj'][item] - update_factor
        else:
            update_cluster_info['ntr_m_Cj'][item] = cluster_info['ntr_m_Cj'][item]
            update_cluster_info['ntr_m_D_Cj'][item] = cluster_info['ntr_m_D_Cj'][item]

    return update_cluster_info


def cluster_info_init(basket, nbaskets, nitems):
    cluster_info = dict()
    cluster_info['nj'] = 1
    cluster_info['ntr_m_Cj'] = [0] * nitems
    cluster_info['ntr_m_D_Cj'] = [0] * nitems

    for item in range(0, nitems):
        if basket[item] == 1:
            cluster_info['ntr_m_Cj'][item] = 1
            cluster_info['ntr_m_D_Cj'][item] = nbaskets - 1
        else:
            cluster_info['ntr_m_Cj'][item] = 0
            cluster_info['ntr_m_D_Cj'][item] = nbaskets

    return cluster_info


class Practical:

    def __init__(self):
        self.clustering = list()
        self.iter_count = 0
        self.cluster_id = 0
        self.max_iter = 100
        self.max_number_of_clusters = 100

    def fit(self, baskets, nbaskets, nitems):

        self.clustering = list()
        self.iter_count = 0
        self.cluster_id = 0

        self.nbaskets = nbaskets
        self.nitems = nitems

        # print datetime.datetime.now(), '_allocation'
        clusters, clusters_info, clusters_quality_dict, clusters_reverse = self._allocation(baskets)

        # print datetime.datetime.now(), '_refinement'
        clusters = self._refinement(baskets, clusters, clusters_info,
                                    clusters_quality_dict, clusters_reverse)

        for cid, cluster in clusters.iteritems():
            self.clustering.append({
                'cluster': cluster,
                'centroid': None
            })
        return self

    def _allocation(self, baskets):

        self.iter_count += 1

        item_baskets = calculate_item_baskets(baskets, self.nbaskets)

        self.w_m_D = calculate_w_m_D(item_baskets, self.nbaskets, self.nitems)

        clusters = defaultdict(dict)
        clusters_info = defaultdict(dict)

        first_element = random.choice(baskets.keys())

        clusters[self.cluster_id][first_element] = baskets[first_element]
        clusters_info[self.cluster_id] = cluster_info_init(baskets[first_element], self.nbaskets, self.nitems)

        self.cluster_id += 1

        clusters_quality_dict = calculate_quality_dict(clusters_info, self.w_m_D, self.nbaskets, self.nitems)

        for bid, basket in baskets.iteritems():
            # print bid, len(clusters)

            if bid == first_element:
                continue

            max_quality = -float('infinity')
            best_cluster_id = None
            best_update = None
            best_cluster_quality = None

            for cid in clusters:

                cur_cluster_info = clusters_info[cid]
                cur_cluster_cid_quality = clusters_quality_dict[cid]
                cluster_info_update = calculate_cluster_info_update(
                        cur_cluster_info, self.w_m_D, self.nbaskets, self.nitems, basket)
                # print cluster_info_update
                cluster_quality_update = calculate_quality_cluster(
                        cluster_info_update, self.w_m_D, self.nbaskets, self.nitems)
                clusters_quality_dict[cid] = cluster_quality_update

                quality_update = calculate_quality(clusters_quality_dict, self.nbaskets)

                clusters_quality_dict[cid] = cur_cluster_cid_quality

                if quality_update > max_quality:
                    max_quality = quality_update
                    best_cluster_id = cid
                    best_update = cluster_info_update
                    best_cluster_quality = cluster_quality_update

            new_cid = self.cluster_id
            self.cluster_id += 1
            cluster_info_new = cluster_info_init(baskets[bid], self.nbaskets, self.nitems)
            quality_new_cluster = calculate_quality_cluster(
                        cluster_info_new, self.w_m_D, self.nbaskets, self.nitems)
            clusters_quality_dict[new_cid] = quality_new_cluster
            new_quality = calculate_quality(clusters_quality_dict, self.nbaskets)

            if len(clusters) < self.max_number_of_clusters and new_quality > max_quality:
                clusters[new_cid][bid] = baskets[bid]
                clusters_info[new_cid] = cluster_info_new
            else:
                clusters[best_cluster_id][bid] = baskets[bid]
                del clusters_quality_dict[new_cid]

                clusters_info[best_cluster_id] = best_update
                clusters_quality_dict[best_cluster_id] = best_cluster_quality

        clusters_reverse = dict()
        for cid, cluster in clusters.iteritems():
            for bid in cluster:
                clusters_reverse[bid] = cid

        return clusters, clusters_info, clusters_quality_dict, clusters_reverse

    def _refinement(self, baskets, clusters, clusters_info, clusters_quality_dict, clusters_reverse):

        while True:
            # print self.iter_count, len(clusters)
            self.iter_count += 1
            moved = False

            cur_quality = calculate_quality(clusters_quality_dict, self.nbaskets)

            for bid, basket in baskets.iteritems():

                # print bid, len(clusters)

                original_cid = clusters_reverse[bid]
                original_cluster_info = clusters_info[original_cid]
                original_cluster_quality = clusters_quality_dict[original_cid]
                cluster_info_update_rem = calculate_cluster_info_update(original_cluster_info,
                                                                        self.w_m_D, self.nbaskets, self.nitems, basket,
                                                                        update_factor=-1)
                cluster_quality_update_rem = None
                if cluster_info_update_rem is not None:
                    cluster_quality_update_rem = calculate_quality_cluster(
                            cluster_info_update_rem, self.w_m_D, self.nbaskets, self.nitems)
                    clusters_quality_dict[original_cid] = cluster_quality_update_rem
                else:
                    del clusters_quality_dict[original_cid]

                max_quality = cur_quality
                best_cluster_id = None
                best_update = None
                best_cluster_quality = None

                for cid in clusters:
                    if cid == original_cid:
                        continue

                    cur_cluster_info = clusters_info[cid]
                    cur_cluster_quality = clusters_quality_dict[cid]
                    cluster_info_update_add = calculate_cluster_info_update(cur_cluster_info,
                                                                        self.w_m_D, self.nbaskets, self.nitems, basket)
                    cluster_quality_update_add = calculate_quality_cluster(
                            cluster_info_update_add, self.w_m_D, self.nbaskets, self.nitems)
                    clusters_quality_dict[cid] = cluster_quality_update_add

                    quality_update = calculate_quality(clusters_quality_dict, self.nbaskets)

                    clusters_quality_dict[cid] = cur_cluster_quality

                    if quality_update > max_quality:
                        max_quality = quality_update
                        best_cluster_id = cid
                        best_update = cluster_info_update_add
                        best_cluster_quality = cluster_quality_update_add

                new_cid = self.cluster_id
                self.cluster_id += 1

                cluster_info_new = cluster_info_init(baskets[bid], self.nbaskets, self.nitems)
                quality_new_cluster = calculate_quality_cluster(
                            cluster_info_new, self.w_m_D, self.nbaskets, self.nitems)
                clusters_quality_dict[new_cid] = quality_new_cluster
                new_quality = calculate_quality(clusters_quality_dict, self.nbaskets)

                if len(clusters) < self.max_number_of_clusters and new_quality > max_quality:

                    moved = True

                    clusters[new_cid][bid] = baskets[bid]
                    clusters_info[new_cid] = cluster_info_new
                    clusters_quality_dict[new_cid] = quality_new_cluster

                    clusters_info[original_cid] = cluster_info_update_rem
                    clusters_quality_dict[original_cid] = cluster_quality_update_rem

                    del clusters[original_cid][bid]

                    if len(clusters[original_cid]) == 0:
                        del clusters[original_cid]
                        del clusters_info[original_cid]
                        del clusters_quality_dict[original_cid]

                    clusters_reverse[bid] = new_cid

                elif max_quality > cur_quality:

                    moved = True

                    clusters[best_cluster_id][bid] = baskets[bid]
                    clusters_info[best_cluster_id] = best_update
                    clusters_quality_dict[best_cluster_id] = best_cluster_quality

                    clusters_info[original_cid] = cluster_info_update_rem
                    clusters_quality_dict[original_cid] = cluster_quality_update_rem

                    del clusters_quality_dict[new_cid]
                    del clusters[original_cid][bid]

                    if len(clusters[original_cid]) == 0:
                        del clusters[original_cid]
                        del clusters_info[original_cid]
                        del clusters_quality_dict[original_cid]

                    clusters_reverse[bid] = best_cluster_id

                else:
                    moved = False

                    del clusters_quality_dict[new_cid]

                    clusters_info[original_cid] = original_cluster_info
                    clusters_quality_dict[original_cid] = original_cluster_quality

            if self.iter_count >= self.max_iter:
                break

            if not moved:
                break

        return clusters
