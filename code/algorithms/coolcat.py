import copy
import random
import numpy as np
from collections import defaultdict

__author__ = 'Riccardo Guidotti'


def entropy_cluster(cluster_info, N):
    n = cluster_info['n']
    if n == 0:
        return 0.0
    item_freq = cluster_info['item_freq']
    entropy_cluster_val = 0
    for item in range(0, len(item_freq)):
        p = 1.0 * item_freq[item] / n
        if p == 0 or p == 1.0:
            entropy_item = 0.0
        else:
            entropy_item = - ((p * np.log2(p)) + ((1-p) * np.log2(1-p)))
        entropy_cluster_val += entropy_item

    return entropy_cluster_val / N


def entropy(clusters_entropy):
    return np.sum(clusters_entropy.values())


def cluster_info_init(basket, nitems):
    cluster_info = dict()

    cluster_info['n'] = 1
    cluster_info['item_freq'] = [0] * nitems

    for item in range(0, len(basket)):
        if basket[item]:
            cluster_info['item_freq'][item] += 1

    return cluster_info


def add_basket(cluster_info, basket):

    cluster_info['n'] += 1

    for item in range(0, len(basket)):
        if basket[item]:
            cluster_info['item_freq'][item] += 1

    return cluster_info


def remove_basket(cluster_info, basket):

    cluster_info['n'] -= 1

    if cluster_info['n'] == 0:
        return None

    for item in range(0, len(basket)):
        if basket[item]:
            cluster_info['item_freq'][item] -= 1

    return cluster_info


class Coolcat:

    def __init__(self):
        self.clustering = list()
        self.iter_count = 0
        self.cluster_id = 0

    def fit(self, baskets, nbaskets, nitems, k, random_sample=float('infinity')):

        self.clustering = list()
        self.iter_count = 0
        self.cluster_id = 0

        self.nbaskets = nbaskets
        self.nitems = nitems
        self.k = k

        self.random_sample = random_sample

        if self.random_sample > self.nbaskets:
            self.random_sample = self.nbaskets

        # print 'initialization'
        clusters, clusters_info, clusters_entropy, seeds_set = self._initialization(baskets)

        # print 'allocation'
        clusters, clusters_info, clusters_entropy, clusters_reverse = self._allocation(
                baskets, clusters, clusters_info, clusters_entropy, seeds_set)

        # print 'refinement'
        clusters = self._refinement(baskets, clusters, clusters_info, clusters_entropy, clusters_reverse)

        # print clusters

        for cluster in clusters.values():
            res = dict()
            for bid in cluster:
                res[bid] = baskets[bid]
            self.clustering.append({
                'cluster': res,
                'centroid': None
            })

        return self

    def _initialization(self, baskets):

        self.iter_count += 1

        baskets_sample_keys = random.sample(baskets.keys(), self.random_sample)

        b1b2_dist = defaultdict(dict)

        max_dist = -float('infinity')
        best_couple = None

        for b1 in baskets_sample_keys:
            cluster_info = cluster_info_init(baskets[b1], self.nitems)
            for b2 in baskets_sample_keys:
                if b1 == b2:
                    continue

                cluster_info = add_basket(cluster_info, baskets[b2])
                dist = entropy_cluster(cluster_info, self.nbaskets)
                b1b2_dist[b1][b2] = dist

                if dist > max_dist:
                    max_dist = dist
                    best_couple = (b1, b2)

                cluster_info = remove_basket(cluster_info, baskets[b2])

        clusters = defaultdict(list)
        clusters_info = dict()
        clusters_entropy = dict()

        cid = self.cluster_id
        clusters[cid].append(best_couple[0])
        clusters_info[cid] = cluster_info_init(baskets[best_couple[0]], self.nitems)
        clusters_entropy[cid] = 0.0
        self.cluster_id += 1

        cid = self.cluster_id
        clusters[cid].append(best_couple[1])
        clusters_info[cid] = cluster_info_init(baskets[best_couple[1]], self.nitems)
        clusters_entropy[cid] = 0.0
        self.cluster_id += 1

        seeds_set = [best_couple[0], best_couple[1]]
        seeds_set = set(seeds_set)

        while len(seeds_set) < self.k:

            max_dist = -float('infinity')
            best_candidate = None
            for s in seeds_set:
                min_dist = float('infinity')
                best_candidate_b = None

                for b in b1b2_dist[s]:
                    if b not in seeds_set:

                        dist = b1b2_dist[b][s]

                        if dist > 0 and dist < min_dist:
                            min_dist = dist
                            best_candidate_b = b

                if min_dist > max_dist:
                    max_dist = min_dist
                    best_candidate = best_candidate_b

            seeds_set.add(best_candidate)
            cid = self.cluster_id
            clusters[cid].append(best_candidate)
            clusters_info[cid] = cluster_info_init(baskets[best_candidate], self.nitems)
            clusters_entropy[cid] = 0.0
            self.cluster_id += 1

        return clusters, clusters_info, clusters_entropy, seeds_set

    def _allocation(self, baskets, clusters, clusters_info, clusters_entropy, seeds_set):

        self.iter_count += 1

        clusters_reverse = dict()

        for b in baskets:
            if b not in seeds_set:
                min_entropy = float('infinity')
                best_cid = None
                best_cluster_info = None
                for cid in clusters:
                    cluster_info = copy.deepcopy(clusters_info[cid])
                    cluster_info = add_basket(cluster_info, baskets[b])
                    new_entropy_cluster = entropy_cluster(cluster_info, self.nbaskets)
                    if new_entropy_cluster < min_entropy:
                        min_entropy = new_entropy_cluster
                        best_cid = cid
                        best_cluster_info = cluster_info

                clusters[best_cid].append(b)
                clusters_info[best_cid] = best_cluster_info
                clusters_entropy[best_cid] = min_entropy
                clusters_reverse[b] = best_cid
            else:
                for cid in clusters:
                    if b in clusters[cid]:
                        clusters_reverse[b] = cid
                        break

        return clusters, clusters_info, clusters_entropy, clusters_reverse

    def _refinement(self, baskets, clusters, clusters_info, clusters_entropy, clusters_reverse):

        while True:
            moved = False
            self.iter_count += 1

            for b in baskets:
                current_entropy = entropy(clusters_entropy)

                original_cid = clusters_reverse[b]

                original_cluster_info = copy.deepcopy(clusters_info[original_cid])
                original_cluster_entropy = clusters_entropy[original_cid]

                clusters_info_tmp = remove_basket(clusters_info[original_cid], baskets[b])
                if clusters_info_tmp is None:
                    clusters_info[original_cid] = original_cluster_info
                    clusters_entropy[original_cid] = original_cluster_entropy
                    continue

                clusters_entropy[original_cid] = clusters_info_tmp
                clusters_entropy[original_cid] = entropy_cluster(clusters_info[original_cid], self.nbaskets)

                for cid in clusters:
                    if cid == original_cid:
                        continue

                    cid_cluster_info = copy.deepcopy(clusters_info[cid])
                    cid_cluster_entropy = clusters_entropy[cid]

                    clusters_info[cid] = add_basket(clusters_info[cid], baskets[b])
                    clusters_entropy[cid] = entropy_cluster(clusters_info[cid], self.nbaskets)

                    new_entropy = entropy(clusters_entropy)

                    if new_entropy < current_entropy:
                        moved = True
                    else:
                        clusters_info[cid] = cid_cluster_info
                        clusters_entropy[cid] = cid_cluster_entropy

                if not moved:
                    clusters_info[original_cid] = original_cluster_info
                    clusters_entropy[original_cid] = original_cluster_entropy

            if not moved:
                break

        return clusters
