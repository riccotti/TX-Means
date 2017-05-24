import random

from heapq import *

from util import *

__author__ = 'Riccardo Guidotti'


def calculate_distances(baskets, centroids):
    res = dict()
    distances = defaultdict(dict)

    for b in baskets:
        res[b] = list()
        for c in centroids:
            dist = jaccard_bitarray(baskets[b], centroids[c])
            heappush(res[b], (dist, c))
            distances[b][c] = dist

    return res, distances


def fun(baskets, m):
    res = 0
    for b in baskets:
        res += jaccard_bitarray(m, baskets[b])
    return 1.0 * res / len(baskets)


def rep(baskets, item_baskets, freq, nbaskets, nitems, min_item_freq=2):

    max_freq = max(freq.values())
    min_freq = min(freq.values())
    m = nitems * bitarray('0')

    if max_freq < min_item_freq:
        return m

    if max_freq == min_freq:
        for e in freq:
            m[e] = 1
        return m

    sorted_freq = sorted(freq, key=freq.get, reverse=True)

    len_m = 0
    for e in sorted_freq:
        if freq[e] < min_item_freq:
            break

        if freq[e] < max_freq:
            max_freq = freq[e]
            break
        m[e] = 1
        len_m += 1

    fun_val1 = fun(baskets, m)

    iterations = 0

    while True and max_freq >= min_freq and not m.all():
        iterations += 1

        m2 = bitarray(m)
        len_m2 = len_m
        nbr2 = nbaskets * bitarray('0')
        for e in sorted_freq[len_m:len(sorted_freq)]:

            if freq[e] < min_item_freq:
                break

            if freq[e] < max_freq:
                max_freq = freq[e]
                break
            m2[e] = 1
            nbr2 |= item_baskets[e]
            len_m2 += 1

        fun_val2 = fun(baskets, m2)

        if fun_val1 <= fun_val2:
            return m
        else:
            m = bitarray(m2)
            nbr = nbr2
            fun_val1 = fun_val2
            len_m = len_m2

    return m


def cluster_baskets(baskets, centroids, distances):

    clusters = dict()

    for c in centroids:
        clusters[c] = list()

    cluster_centroid = dict()

    for b in baskets:
        min_val = heappop(distances[b])
        c = min_val[1]
        clusters[c].append(b)
        cluster_centroid[b] = c

    return clusters, cluster_centroid


def equals_centroids(centroids, centroids_new):

    for e1 in centroids.values():
        if e1 not in centroids_new.values():
            return False

    for e2 in centroids_new.values():
        if e2 not in centroids.values():
            return False

    return True


def calculate_item_baskets_freq(baskets, nbaskets):
    freq = defaultdict(int)
    item_baskets = defaultdict(lambda: nbaskets * bitarray('0'))
    for b in baskets:
        for item in range(0, len(baskets[b])):
            if baskets[b][item]:
                item_baskets[item][b] = 1
                freq[item] += 1
    return item_baskets, freq


class TKMeans:

    def __init__(self, ):

        self.clustering = list()
        self.iter_count = 0

    def fit(self, baskets, nbaskets, nitems, k, min_item_freq=2, niter=100, max_iter=100):

        self.clustering = list()
        self.iter_count = 0

        self.nbaskets = nbaskets
        self.nitems = nitems
        self.k = k
        self.min_item_freq = min_item_freq
        self.niter = niter
        self.max_iter = max_iter

        max_dist = float('infinity')
        best_res = None
        best_iter_count = None
        for i in range(0, self.niter):

            clustering_res, distances, iter_count = self._run(baskets)
            tot_dist = 0
            den = 0

            for cluster_dict in clustering_res:
                c_id = cluster_dict['centroid_id']

                if len(cluster_dict['cluster']) > 0:
                    for b in cluster_dict['cluster']:
                        tot_dist += distances[b][c_id]
                        den += 1
                    tot_dist /= den

            if tot_dist < max_dist:
                max_dist = tot_dist
                best_res = clustering_res
                best_iter_count = iter_count

        for cluster in best_res:

            self.clustering.append({
                        'cluster': cluster['cluster'],
                        'centroid': cluster['centroid']
                    })

        self.iter_count = best_iter_count

        return self

    def _run(self, baskets):

        selected_seeds = random.sample(baskets.keys(), self.k)

        centroids = dict()
        for s in selected_seeds:
            centroids[s] = baskets[s]

        clusters = None
        distances = None
        cluster_centroid = None
        iter_count = 0
        while True:
            if iter_count >= self.max_iter:
                break
            # print 'iter', iter
            iter_count += 1

            dist, distances = calculate_distances(baskets, centroids)
            clusters, cluster_centroid = cluster_baskets(baskets, centroids, dist)

            new_centroids = dict()
            for c, i in zip(clusters, range(1, len(clusters)+1)):
                subbaskets = dict()
                for b in clusters[c]:
                    subbaskets[b] = baskets[b]

                if len(subbaskets) > 0:
                    item_baskets, freq = calculate_item_baskets_freq(subbaskets, self.nbaskets)
                    new_centroids[-i*1000] = rep(subbaskets, item_baskets, freq, self.nbaskets,
                                                 self.nitems, self.min_item_freq)
                else:
                    new_centroids[-i*1000] = self.nitems * bitarray('0')

            if equals_centroids(centroids, new_centroids):
                break

            centroids = new_centroids

        clustering_res = list()

        for cluster in clusters.values():
            res = dict()
            for c in cluster:
                res[c] = baskets[c]

            if len(res) > 0:
                clustering_res.append({
                    'cluster': res,
                    'centroid': centroids[cluster_centroid[res.keys()[0]]],
                    'centroid_id': cluster_centroid[res.keys()[0]]
                })
            else:
                clustering_res.append({
                    'cluster': dict(),
                    'centroid': self.nitems * bitarray('0'),
                    'centroid_id': None
                })

        return clustering_res, distances, iter_count
