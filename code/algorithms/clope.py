import copy
import random

from bitarray import bitarray
from collections import defaultdict

__author__ = 'Riccardo Guidotti'


def delta_profit(cluster_info, basket, r):

    if cluster_info is None:
        N_new = 1
        S_new = basket.count()
        W_new = basket.count()
        return 1.0 * S_new * N_new / W_new ** r

    N_new = cluster_info['N'] + 1
    S_new = cluster_info['S'] + basket.count()
    W_new = (cluster_info['set'] | basket).count()
    P_new = 1.0 * S_new * N_new / W_new ** r
    P_old = 1.0 * cluster_info['S'] * cluster_info['N'] / cluster_info['W'] ** r
    return P_new - P_old


def delta_profit_minus(cluster_info, basket, r, cluster):

    N_new = cluster_info['N'] - 1
    S_new = cluster_info['S'] - basket.count()
    tmp = bitarray('0') * len(basket)
    for bbasket in cluster.values():
        tmp |= bbasket
    W_new = tmp.count()

    P_new = 1.0 * S_new * N_new / W_new ** r
    P_old = 1.0 * cluster_info['S'] * cluster_info['N'] / cluster_info['W'] ** r
    return P_new - P_old


def cluster_profit(cluster_info, r):
    return (1.0 * cluster_info['S'] / cluster_info['W'] ** r) * cluster_info['N']


def cluster_info_init(basket, nbaskets, nitems):
    cluster_info = dict()

    cluster_info['N'] = 1
    cluster_info['S'] = basket.count()
    cluster_info['set'] = copy.deepcopy(basket)
    cluster_info['W'] = cluster_info['set'].count()

    return cluster_info


class Clope:

    def __init__(self):
        self.clustering = list()
        self.iter_count = 0
        self.cluster_id = 0
        self.max_number_of_clusters = 100

    def fit(self, baskets, nbaskets, nitems, r):

        self.clustering = list()
        self.iter_count = 0
        self.cluster_id = 0

        self.nbaskets = nbaskets
        self.nitems = nitems
        self.r = r

        # print 'allocation'
        clusters, clusters_reverse, clusters_info = self._allocation(baskets)

        # for cid, cluster in clusters.iteritems():
        #     print cid, cluster, clusters_info[cid]

        # print 'refinement'
        clusters = self._refinement(baskets, clusters, clusters_reverse, clusters_info)

        for cid, cluster in clusters.iteritems():
            self.clustering.append({
                'cluster': cluster,
                'centroid': None
            })
        return self

    def _allocation(self, baskets):

        self.iter_count += 1

        clusters = defaultdict(dict)
        clusters_info = defaultdict(dict)
        clusters_reverse = dict()

        first_element = random.choice(baskets.keys())

        clusters[self.cluster_id][first_element] = baskets[first_element]
        clusters_info[self.cluster_id] = cluster_info_init(baskets[first_element], self.nbaskets, self.nitems)
        clusters_reverse[first_element] = self.cluster_id

        current_profit_num = cluster_profit(clusters_info[self.cluster_id], self.r)
        den = 1

        self.cluster_id += 1

        for bid, basket in baskets.iteritems():

            if bid == first_element:
                continue

            max_profit = -float('infinity')
            best_cluster_id = None
            for cid in clusters:

                d_profit = delta_profit(clusters_info[cid], basket, self.r)
                profit_in_cluster = (current_profit_num + d_profit) / (den + 1)

                if profit_in_cluster > max_profit:
                    max_profit = profit_in_cluster
                    best_cluster_id = cid

            d_profit = delta_profit(None, basket, self.r)
            profit_new_cluster = (current_profit_num + d_profit) / (den + 1)

            if len(clusters) < self.max_number_of_clusters and profit_new_cluster > max_profit:
                new_cid = self.cluster_id
                self.cluster_id += 1
                clusters[new_cid][bid] = baskets[bid]
                clusters_info[new_cid] = cluster_info_init(baskets[bid], self.nbaskets, self.nitems)
                clusters_reverse[bid] = new_cid
                current_profit_num += delta_profit(None, basket, self.r)
            else:
                # print current_profit_num, delta_profit(clusters_info[best_cluster_id], basket, self.r), 'delta'
                current_profit_num += delta_profit(clusters_info[best_cluster_id], basket, self.r)
                clusters[best_cluster_id][bid] = baskets[bid]
                clusters_info[best_cluster_id]['N'] += 1
                clusters_info[best_cluster_id]['S'] += basket.count()
                clusters_info[best_cluster_id]['set'] |= basket
                clusters_info[best_cluster_id]['W'] = clusters_info[best_cluster_id]['set'].count()
                clusters_reverse[bid] = best_cluster_id
            den += 1

        # print current_profit_num, '<<<<<<'
        # for cid in clusters_info:
        #     print cid, cluster_profit(clusters_info[cid], self.r)
        return clusters, clusters_reverse, clusters_info

    def _refinement(self, baskets, clusters, clusters_reverse, clusters_info):

        current_profit_num = 0
        den = 0
        for cid in clusters_info:
            # print cid, cluster_profit(clusters_info[cid], self.r)
            current_profit_num += cluster_profit(clusters_info[cid], self.r)
            den += clusters_info[cid]['N']

        # print 'INIZIO', current_profit_num

        while True:
            # print ''
            # print self.iter_count, len(clusters)

            self.iter_count += 1
            moved = False

            for bid, basket in baskets.iteritems():

                # print bid, len(clusters)

                current_profit = 1.0 * current_profit_num / den
                # print 'quassu', current_profit_num, 1.0*current_profit_num / den

                original_cid = clusters_reverse[bid]

                best_cid = None
                max_profit = -float('infinity')

                tmp_profit_num = current_profit_num + delta_profit_minus(clusters_info[original_cid], basket, self.r,
                                                                         clusters[original_cid])

                for cid in clusters:
                    if cid == original_cid:
                        continue

                    # print clusters_info[cid]
                    d_profit = delta_profit(clusters_info[cid], basket, self.r)
                    profit_in_cluster = (tmp_profit_num + d_profit) / den

                    # print profit_in_cluster, d_profit, '<-----------'

                    if profit_in_cluster > max_profit:
                        max_profit = profit_in_cluster
                        best_cid = cid

                d_profit = delta_profit(None, basket, self.r)
                profit_new_cluster = (tmp_profit_num + d_profit) / den

                # print bid, profit_new_cluster, max_profit, current_profit

                if profit_new_cluster > max_profit and profit_new_cluster > current_profit:

                    d_profit = delta_profit(None, basket, self.r)
                    current_profit_num = tmp_profit_num + d_profit

                    new_cid = self.cluster_id
                    self.cluster_id += 1
                    clusters[new_cid][bid] = baskets[bid]
                    clusters_info[new_cid] = cluster_info_init(baskets[bid], self.nbaskets, self.nitems)

                    clusters_reverse[bid] = new_cid

                    # print 'new, remove ', bid, ' from ', original_cid, ' put in ', new_cid

                    del clusters[original_cid][bid]
                    if len(clusters[original_cid]) == 0:
                        del clusters[original_cid]
                        del clusters_info[original_cid]
                    else:
                        clusters_info[original_cid]['N'] -= 1
                        clusters_info[original_cid]['S'] -= basket.count()
                        clusters_info[original_cid]['set'] = bitarray('0') * self.nitems
                        for bbasket in clusters[original_cid].values():
                            clusters_info[original_cid]['set'] |= bbasket
                        clusters_info[original_cid]['W'] = clusters_info[original_cid]['set'].count()

                    moved = True
                    clusters_reverse[bid] = new_cid

                elif max_profit > current_profit:

                    # print 'move, remove ', bid, ' from ', original_cid, ' put in ', best_cid

                    d_profit = delta_profit(clusters_info[best_cid], basket, self.r)
                    current_profit_num = tmp_profit_num + d_profit

                    clusters[best_cid][bid] = baskets[bid]
                    clusters_info[best_cid]['N'] += 1
                    clusters_info[best_cid]['S'] += basket.count()
                    clusters_info[best_cid]['set'] |= basket
                    clusters_info[best_cid]['W'] = clusters_info[best_cid]['set'].count()

                    del clusters[original_cid][bid]
                    if len(clusters[original_cid]) == 0:
                        del clusters[original_cid]
                        del clusters_info[original_cid]
                    else:
                        clusters_info[original_cid]['N'] -= 1
                        clusters_info[original_cid]['S'] -= basket.count()
                        clusters_info[original_cid]['set'] = bitarray('0') * self.nitems
                        for bbasket in clusters[original_cid].values():
                            clusters_info[original_cid]['set'] |= bbasket
                        clusters_info[original_cid]['W'] = clusters_info[original_cid]['set'].count()

                    moved = True
                    clusters_reverse[bid] = best_cid

            if not moved:
                break

            # if self.iter_count == 10:
            #     break

        return clusters


# def main():
#     print 'Test'
#     baskets_list = [
#     [5, 4, 1, 3, 6, 2] ,
#     [1, 0, 5, 3, 6, 4] ,
#     [10, 11, 7, 9, 8, 6] ,
#     [48, 6, 10, 7, 11, 8] ,
#     [18, 15, 14, 17, 12, 16] ,
#     [17, 12, 13, 16, 14, 18] ,
#     [49, 23, 20, 19, 21, 18] ,
#     [23, 19, 22, 20, 18, 49],
#     [26, 24, 33, 25, 34, 28, 31, 35, 30, 32, 29, 27] ,
#     [27, 26, 24, 29, 34, 33, 35, 30, 31, 28, 32, 25] ,
#     [30, 25, 33, 35, 27, 26, 34, 32, 24, 31, 29, 28] ,
#     [24, 31, 33, 34, 26, 35, 32, 25, 29, 27, 30, 28] ,
#     [46, 39, 36, 41, 43, 45, 44, 40, 47, 37, 42, 38] ,
#     [47, 42, 38, 46, 40, 37, 36, 45, 43, 41, 39, 44] ,
#     [43, 39, 46, 47, 45, 36, 40, 44, 42, 37, 38, 41] ,
#     [43, 45, 36, 46, 39, 47, 37, 42, 38, 40, 41, 44]
#     ]
#
#     # baskets_list = [
#     # [1,2,3,4,5,6] , [1,2,3,4,5,6],
#     # [2,3,4,5,6,7], [2,3,4,5,6,7]
#     # ]
#
#     # baskets_list = [
#     #     [1, 2, 3], [1, 2, 3, 4],
#     #     [2, 3, 4, 5], [3, 4, 5]
#     # ]
#
#     # baskets_list = [
#     #     [1, 2, 3], [3, 4, 5],
#     #     [6, 7, 8], [8, 9, 10]
#     # ]
#
#
#     baskets_list, map_newitem_item, map_item_newitem = remap_items(baskets_list)
#     baskets_list = basket_list_to_bitarray(baskets_list, len(map_newitem_item))
#     practical = Clope(r=2)
#     practical.fit(baskets_list)
#
#     print '------------------'
#     for c in practical.clustering:
#         print c['cluster'].keys()
#
# if __name__ == "__main__":
#     main()