import copy
from heapq import *
from bitarray import bitarray
from collections import defaultdict

__author__ = 'Riccardo Guidotti'


def calculate_item_baskets(baskets, nbaskets, nitems):
    item_baskets = defaultdict(lambda: nbaskets * bitarray('0'))
    N_a_dict = [0] * nitems
    for b in baskets:
        for item in range(0, len(baskets[b])):
            if baskets[b][item]:
                item_baskets[item][b] = 1
                N_a_dict[item] += 1

    return item_baskets, N_a_dict


def calc_clusters_info(item_baskets, nbaskets, nitems):
    cluster_info = dict()
    cluster_info['n'] = nbaskets
    cluster_info['n_a'] = [0] * nitems

    for item in item_baskets:
        cluster_info['n_a'][item] = item_baskets[item].count()

    return cluster_info


def quality_cluster(cluster_info, N_a_dict, N, nitems):
    quality_val = 0
    n = cluster_info['n']
    for item in range(0, nitems):
        n_a = cluster_info['n_a'][item]
        N_a = N_a_dict[item]
        quality_val += ((1.0 * n_a / n)**2 - (1.0 * N_a / N)**2)

    return quality_val * n / N


def quality(clusters_quality, N):
    return 1.0 * sum(clusters_quality.values())


def update_cluster_info(basket, cluster_info, N_a_dict, sign, N, nitems):

    cluster_info_new = dict()
    cluster_info_new['n'] = cluster_info['n'] + 1.0 * sign
    cluster_info_new['n_a'] = [0] * nitems

    quality_val_new = 0
    n = cluster_info_new['n']

    for item in range(0, nitems):
        if basket[item]:
            cluster_info_new['n_a'][item] = cluster_info['n_a'][item] + 1.0 * sign
        else:
            cluster_info_new['n_a'][item] = cluster_info['n_a'][item]

        n_a = cluster_info_new['n_a'][item]
        N_a = N_a_dict[item]
        if n > 0 and n_a > 0:
            quality_val_new += ((1.0 * n_a / n)**2 - (1.0 * N_a / N)**2)

    return cluster_info_new, quality_val_new * n / N


class Atdc:

    def __init__(self):
        self.clustering = list()
        self.iter_count = 0
        self.cluster_id = 0
        self.max_iter = 100
        self.max_number_of_clusters = 100

    def fit(self, baskets, nbaskets, nitems, min_cluster_size=2):

        self.clustering = list()
        self.iter_count = 0
        self.cluster_id = 0

        self.nbaskets = nbaskets
        self.nitems = nitems
        self.min_cluster_size = min_cluster_size

        clusters = dict()
        clusters_quality = dict()
        clusters_info = dict()
        clusters_quality_list = list()
        basket_cluster = dict()
        self.item_baskets_all, self.N_a_dict = calculate_item_baskets(baskets, self.nbaskets, self.nitems)

        clusters[self.cluster_id] = bitarray('1') * self.nbaskets
        clusters_info[self.cluster_id] = calc_clusters_info(self.item_baskets_all, self.nbaskets, self.nitems)
        clusters_quality[self.cluster_id] = quality_cluster(
                clusters_info[self.cluster_id], self.N_a_dict, self.nbaskets, self.nitems)
        for bid in baskets:
            basket_cluster[bid] = self.cluster_id

        heappush(clusters_quality_list, (clusters_quality[self.cluster_id], self.cluster_id))

        cur_quality = quality(clusters_quality, self.nbaskets)

        self.cluster_id += 1

        self.iter_count += 1

        while True:
            # print 'iter', self.iter_count, len(clusters)
            self.iter_count += 1

            moved = False

            new_cid = self.cluster_id

            self.cluster_id += 1

            clusters[new_cid] = bitarray('0') * self.nbaskets
            clusters_info[new_cid] = calc_clusters_info(dict(), 0, self.nitems)
            clusters_quality[new_cid] = 0.0

            while True:
                # print clusters_quality_list
                # print clusters
                min_val = heappop(clusters_quality_list)
                cur_cid = min_val[1]
                # print 'cid to be split', cur_cid, new_cid

                clusters_new = copy.deepcopy(clusters)
                clusters_info_new = copy.deepcopy(clusters_info)
                clusters_quality_new = copy.deepcopy(clusters_quality)
                basket_cluster_new = copy.deepcopy(basket_cluster)

                # print 'partitioning'
                clusters_new, clusters_info_new, clusters_quality_new, basket_cluster_new = self._partition_cluster(
                        baskets, clusters_new, clusters_info_new, clusters_quality_new, basket_cluster_new,
                        cur_cid, new_cid)

                new_quality = quality(clusters_quality_new, self.nbaskets)

                if len(clusters) < self.max_number_of_clusters and new_quality > cur_quality and \
                                clusters_info_new[new_cid]['n'] >= self.min_cluster_size and \
                                clusters_info_new[cur_cid]['n'] >= self.min_cluster_size:

                    moved = True

                    stabilize = True
                    if clusters_info_new[new_cid]['n'] == 0:
                        stabilize = False
                        del clusters_new[new_cid]
                        del clusters_info_new[new_cid]
                        del clusters_quality_new[new_cid]

                    if clusters_info_new[cur_cid]['n'] == 0:
                        stabilize = False
                        del clusters_new[cur_cid]
                        del clusters_info_new[cur_cid]
                        del clusters_quality_new[cur_cid]

                    if stabilize:
                        # print 'stabilize'
                        clusters_new, clusters_info_new, clusters_quality_new, basket_cluster_new = self._stabilize_clusters(
                            baskets, clusters_new, clusters_info_new, clusters_quality_new, basket_cluster_new)

                    clusters = clusters_new
                    clusters_info = clusters_info_new
                    clusters_quality = clusters_quality_new
                    basket_cluster = basket_cluster_new

                    if new_cid in clusters:
                        q_new = clusters_quality[new_cid]
                        heappush(clusters_quality_list, (q_new, new_cid))

                    if cur_cid in clusters:
                        q_cid = clusters_quality[cur_cid]
                        heappush(clusters_quality_list, (q_cid, cur_cid))

                    cur_quality = quality(clusters_quality, self.nbaskets)

                    break

                if len(clusters_quality_list) == 0:
                    del clusters[new_cid]
                    del clusters_info[new_cid]
                    del clusters_quality[new_cid]
                    break

            if self.iter_count >= self.max_iter:
                break

            if not moved:
                break

        for cluster in clusters.values():
            res = dict()
            for bid in range(0, len(cluster)):
                if cluster[bid]:
                    res[bid] = baskets[bid]
            self.clustering.append({
                'cluster': res,
                'centroid': None
            })

        return self

    def _partition_cluster(self, baskets, clusters, clusters_info, clusters_quality, basket_cluster, c0, c1):

        list_of_bids = clusters[c0] | clusters[c1]
        while True:
            moved = False

            for bid in range(0, self.nbaskets):

                if not list_of_bids[bid]:
                    continue

                if clusters[c0][bid]:
                    c_u = c0
                    c_v = c1
                else:
                    c_u = c1
                    c_v = c0

                q_u = clusters_quality[c_u]
                q_v = clusters_quality[c_v]
                Qi = q_u + q_v

                c_u_cluster_info_min, q_min_bid = update_cluster_info(
                        baskets[bid], clusters_info[c_u], self.N_a_dict, -1, self.nbaskets, self.nitems)

                c_v_cluster_info_plus, q_plus_bid = update_cluster_info(
                        baskets[bid], clusters_info[c_v], self.N_a_dict, 1, self.nbaskets, self.nitems)

                Qs = q_min_bid + q_plus_bid

                if Qs > Qi:
                    moved = True
                    # print 'qui'
                    clusters[c_v][bid] = 1
                    clusters[c_u][bid] = 0

                    clusters_info[c_u] = c_u_cluster_info_min
                    clusters_info[c_v] = c_v_cluster_info_plus

                    clusters_quality[c_u] = q_min_bid
                    clusters_quality[c_v] = q_plus_bid

                    basket_cluster[bid] = c_v

            if not moved:
                break

        return clusters, clusters_info, clusters_quality, basket_cluster

    def _stabilize_clusters(self, baskets, clusters, clusters_info, clusters_quality, basket_cluster):

        while True:
            moved = False
            cur_quality = quality(clusters_quality, self.nbaskets)
            for bid in baskets:

                c_pivot = basket_cluster[bid]

                cid_list = clusters.keys()
                for cid in cid_list:

                    if cid == c_pivot or cid is None or c_pivot is None:
                        continue

                    c_pivot_cluster_info_min, q_min_bid = update_cluster_info(
                        baskets[bid], clusters_info[c_pivot], self.N_a_dict, -1, self.nbaskets, self.nitems)

                    c_cid_cluster_info_plus, q_plus_bid = update_cluster_info(
                            baskets[bid], clusters_info[cid], self.N_a_dict, 1, self.nbaskets, self.nitems)

                    q_pivot_old = clusters_quality[c_pivot]
                    q_cid_old = clusters_quality[cid]

                    clusters_quality[c_pivot] = q_min_bid
                    clusters_quality[cid] = q_plus_bid

                    new_quality = quality(clusters_quality, self.nbaskets)

                    if new_quality > cur_quality:

                        clusters[cid][bid] = 1
                        clusters_info[cid] = c_cid_cluster_info_plus
                        basket_cluster[bid] = cid

                        if c_pivot_cluster_info_min['n'] == 0:

                            del clusters[c_pivot]
                            del clusters_info[c_pivot]
                            del clusters_quality[c_pivot]

                        clusters[c_pivot][bid] = 0
                        clusters_info[c_pivot] = c_pivot_cluster_info_min
                        c_pivot = cid
                        cur_quality = new_quality

                    else:
                        clusters_quality[c_pivot] = q_pivot_old
                        clusters_quality[cid] = q_cid_old

            if not moved:
                break

        return clusters, clusters_info, clusters_quality, basket_cluster


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
#
#     # baskets_list = [
#     # [1,2,3] ,
#     # [1,2,3] ,
#     # [3,4,5] ,
#     # [3,4,5] ,
#     # [14,15,16],
#     # [14,15,16],
#     # [16,17,18],
#     # [16,17,18],
#     # [6,7,8,9] ,
#     # [6,7,8,9] ,
#     # [6,7,8,9] ,
#     # [6,7,8,9],
#     # [10,11,12,13] ,
#     # [10,11,12,13] ,
#     # [10,11,12,13] ,
#     # [10,11,12,13],
#     # ]
#
#     # baskets_list = [
#     #     [1,2,3,6,8], [1,4,7,9],
#     #     [2,3,5,7,8], [4,5,6,9],
#     #     [6,8], [6,8,9],
#     #     [7,8], [7,9]
#     # ]
#     # baskets_list = [
#     # [1,2,3,4,5,6] , [1,2,3,4,5,6],
#     # [2,3,4,5,6,7], [2,3,4,5,6,7]
#     # ]
#
#     # baskets_list = [
#     #     [1,2], [1,2],
#     #     [3,4], [3,4]
#     # ]
#     #
#     # baskets_list = list()
#     # for i in range(0, 2):
#     #     baskets_list.append([1,2,3])
#     # for i in range(0, 2):
#     #     baskets_list.append([3,4,5])
#     # for i in range(0, 2):
#     #     baskets_list.append([5,6,7])
#
#     # baskets_list = [
#     #         [1,2,3], [1,2,4], [1,2,5], [1,3,4], [1,3,5], [1,4,5],
#     #         [2,3,4], [2,3,5], [2,4,5], [3,4,5],
#     #
#     #         [1,2,6], [1,2,7], [1,6,7], [2,6,7]
#     #     ]
#     # baskets_list = [[1,2,3], [1,2,4], [1,2,5], [3,4,5], [4,5,6]]
#
#
#     baskets_list, map_newitem_item, map_item_newitem = remap_items(baskets_list)
#     baskets_list = basket_list_to_bitarray(baskets_list, len(map_newitem_item))
#     nbaskets = len(baskets_list)
#     nitems = count_items(baskets_list)
#     atdc = Atdc()
#     atdc.fit(baskets_list, nbaskets, nitems, min_cluster_size=2)
#
#     print '------------------'
#     for c in atdc.clustering:
#         print c['cluster'].keys()
#
# if __name__ == "__main__":
#     main()