import random
import datetime
import numpy as np

from scipy.special import binom
from collections import defaultdict

from util import *

__author__ = 'Riccardo Guidotti'


def fun(baskets, m):
    res = 0
    for b in baskets:
        res += jaccard_bitarray(m, baskets[b])
    return 1.0 * res / len(baskets)


def calculate_freq_baskets(baskets):
    freq = defaultdict(int)
    for b in baskets:
        for item in range(0, len(baskets[b])):
            if baskets[b][item]:
                freq[item] += 1
    return freq


def calculate_frequencies(item_baskets):
    freq = defaultdict(int)
    for item in item_baskets:
        freq[item] = item_baskets[item].count()
    return freq


def rep(baskets, item_baskets, freq, nbaskets, nitems, min_item_freq=2):

    max_freq = max(freq.values())
    min_freq = min(freq.values())
    m = nitems * bitarray('0')
    nbr = nbaskets * bitarray('0')

    if max_freq < min_item_freq:
        return m, nbr

    if max_freq == min_freq:
        for e in freq:
            m[e] = 1
            nbr |= item_baskets[e]
        return m, nbr

    sorted_freq = sorted(freq, key=freq.get, reverse=True)

    len_m = 0
    for e in sorted_freq:
        if freq[e] < min_item_freq:
            break

        if freq[e] < max_freq:
            max_freq = freq[e]
            break
        m[e] = 1
        nbr |= item_baskets[e]
        len_m += 1

    fun_val1 = fun(baskets, m)

    iterations = 0

    while True and max_freq >= min_freq and not m.all():
        iterations += 1

        m2 = bitarray(m)
        len_m2 = len_m
        nbr2 = nbaskets * bitarray('0') | nbr
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
            return m, nbr
        else:
            m = bitarray(m2)
            nbr = nbr2
            fun_val1 = fun_val2
            len_m = len_m2

    return m, nbr


def calculate_item_baskets(baskets, nbaskets):
    item_baskets = defaultdict(lambda: nbaskets * bitarray('0'))
    for b in baskets:
        for item in range(0, len(baskets[b])):
            if baskets[b][item]:
                item_baskets[item][b] = 1
    return item_baskets


def calculate_neighbors(baskets, item_baskets, nbaskets, nitems):
    neighbors = defaultdict(lambda: nbaskets * bitarray('0'))
    nodes = nbaskets * bitarray('0')
    for b in baskets:
        nodes[b] = 1
        for item in range(0, nitems):
            if baskets[b][item]:
                neighbors[b] |= item_baskets[item]
        neighbors[b][b] = 0

    no_neighbors = defaultdict(lambda: nbaskets * bitarray('0'))
    for b in baskets:
        no_neighbors[b] = ~neighbors[b] & nodes
        no_neighbors[b][b] = 1

    return neighbors, no_neighbors


def calculate_distance(b1, b2, nbr_b1, nbr_b2, use_neighbors):
    fact1 = jaccard_bitarray(b1, b2)
    dist = fact1
    if use_neighbors:
        fact2 = jaccard_bitarray(nbr_b1, nbr_b2)
        # dist = fact1 * fact2
        # dist = 0.1 * fact1 + 0.9 * fact2
        if fact1 + fact2 > 0:
            dist = 2 * fact1 * fact2 / (fact1 + fact2)
        else:
            dist = 0.0
    return dist


def calculate_distances(baskets, neighbors, m0, m1, nbr_m0, nbr_m1, use_neighbors):
    distances_dict = defaultdict(dict)

    for b in baskets:
        basket = baskets[b]
        if use_neighbors:
            nbr_basket0 = neighbors[0][b]
            nbr_basket1 = neighbors[1][b]
        else:
            nbr_basket0 = None
            nbr_basket1 = None

        distances_dict[-1][b] = calculate_distance(m0, basket, nbr_m0, nbr_basket0, use_neighbors)
        distances_dict[-2][b] = calculate_distance(m1, basket, nbr_m1, nbr_basket1, use_neighbors)

    return distances_dict


def free_params(num_clusters, num_dims):
    return num_clusters * (num_dims + 1)


def cluster_variance(num_points, clusters, distances_dict_values_cl):
    s = 0
    den = float(num_points - len(clusters))
    for distances_dict_values in list(distances_dict_values_cl):
        distances = np.asarray(distances_dict_values)
        s += (distances*distances).sum()

    if 0 != den:
        return s
    else:
        return 0


def loglikelihood(num_points, num_dims, clusters, distances_dict_values_cl):
    ll = 0
    variance = cluster_variance(num_points, clusters, distances_dict_values_cl) or np.nextafter(0, 1)
    # print('var', variance)
    for cluster in clusters:
        fRn = len(cluster)
        t1 = fRn * np.log(fRn)
        t2 = fRn * np.log(num_points)
        t3 = ((fRn * num_dims) / 2.0) * np.log((2.0 * np.pi) * variance)
        t4 = (fRn - 1.0) / 2.0
        ll += t1 - t2 - t3 - t4
    return ll


def bic(clusters, num_points, num_dims, distances_dict_values_cl):

    # num_dims = min(num_points, num_dims)

    log_likelihood = loglikelihood(num_points, num_dims, clusters, list(distances_dict_values_cl))
    num_params = free_params(len(clusters), num_dims)

    return log_likelihood - num_params / 2.0 * np.log(num_points)


def clusters_assignment(baskets, item_baskets, distances_dict, nbaskets, nitems):

    cluster01 = [dict(), dict()]
    freq01 = [defaultdict(int), defaultdict(int)]
    cluster01_assignment = [nbaskets * bitarray('0'), nbaskets * bitarray('0')]
    item_baskets_split = [defaultdict(lambda: nbaskets * bitarray('0')),
                          defaultdict(lambda: nbaskets * bitarray('0'))]

    for b in baskets:
        dist_m0_b = distances_dict[-1][b]
        dist_m1_b = distances_dict[-2][b]
        if dist_m0_b <= dist_m1_b:
            index_cluster = 0
        else:
            index_cluster = 1

        cluster01[index_cluster][b] = baskets[b]
        cluster01_assignment[index_cluster][b] = 1
        for item in range(0, nitems):
            if baskets[b][item]:
                freq01[index_cluster][item] += 1
                item_baskets_split[index_cluster][item][b] = 1

    return cluster01, freq01, cluster01_assignment, item_baskets_split


# def is_homogenous(freq):
#     return np.min(freq.values()) == np.max(freq.values())


def check_convergence_centroids(m0, m1, m0_new, m1_new):
    return (m0_new == m0 and m1_new == m1) or (m0_new == m1 and m1_new == m0)


# def check_convergence_points(cluster01_assignment, cluster01_assignment_new):
#     flag1 = jaccard_bitarray(cluster01_assignment[0], cluster01_assignment_new[0])
#     flag2 = jaccard_bitarray(cluster01_assignment[1], cluster01_assignment_new[1])
#     flag3 = jaccard_bitarray(cluster01_assignment[0], cluster01_assignment_new[1])
#     flag4 = jaccard_bitarray(cluster01_assignment[1], cluster01_assignment_new[0])
#     return (flag1 and flag2) or (flag3 and flag4)


def bisecting_kmeans(baskets, neighbors, item_baskets, distances_dict_init, nbaskets, nitems,
                     m0, m1, use_neighbors, min_item_freq=2, verbose=False):

    if verbose:
        print('\t', datetime.datetime.now(), 'bk initial clusters assignment')

    cluster01, freq01, cluster01_assignment, item_baskets01 = clusters_assignment(
            baskets, item_baskets, distances_dict_init, nbaskets, nitems)

    d0 = distances_dict_init[-1]
    d1 = distances_dict_init[-2]

    if len(freq01[0]) == 0 or len(freq01[1]) == 0:

        if verbose:
            print('\t', datetime.datetime.now(), 'bk the cluster cannot be split')

        if len(cluster01[0]) == 0:
            m0 = m1
            m1 = None
            cluster01[0] = cluster01[1]
            cluster01[1] = None
            d0 = distances_dict_init[-2]
            d1 = None

        if len(cluster01[1]) == 0:
            m0 = m0
            m1 = None
            cluster01[0] = cluster01[0]
            cluster01[1] = None
            d0 = distances_dict_init[-1]
            d1 = None

        bisecting_kmeans_res = {
            'm0': m0,
            'm1': m1,
            'c0': cluster01[0],
            'c1': cluster01[1],
            'd0': d0,
            'd1': None if d1 is None else d1,
            'i': 0
        }

        return bisecting_kmeans_res

    count = 0
    m0_old = None
    m1_old = None

    while True:

        if verbose:
            print('\t', datetime.datetime.now(), 'bk iter', count)

        if count >= 10000:
            break

        count += 1

        neighbors01 = [None, None]

        if use_neighbors:

            if verbose:
                print('\t', datetime.datetime.now(), 'bk calculate neighbors')

            neighbors0, _ = calculate_neighbors(cluster01[0], item_baskets01[0], nbaskets, nitems)
            neighbors1, _ = calculate_neighbors(cluster01[1], item_baskets01[1], nbaskets, nitems)
            neighbors01 = [neighbors0, neighbors1]

        if verbose:
            print('\t', datetime.datetime.now(), 'bk calculate rep m0')

        m0_new, nbr_m0 = rep(cluster01[0], item_baskets01[0], freq01[0], nbaskets, nitems, min_item_freq)

        if verbose:
            print('\t', datetime.datetime.now(), 'bk calculate rep m1')

        m1_new, nbr_m1 = rep(cluster01[1], item_baskets01[1], freq01[1], nbaskets, nitems, min_item_freq)

        if verbose:
            print('\t', datetime.datetime.now(), 'bk check convergence')

        if check_convergence_centroids(m0, m1, m0_new, m1_new):
            break

        if m0_old is not None and check_convergence_centroids(m0_old, m1_old, m0_new, m1_new):
            break

        m0_old = m0
        m1_old = m1

        m0 = m0_new
        m1 = m1_new

        if verbose:
            print('\t', datetime.datetime.now(), 'bk calculate distances with centroids')

        distances_dict = calculate_distances(baskets, neighbors01, m0, m1, nbr_m0, nbr_m1, use_neighbors)

        if verbose:
            print('\t', datetime.datetime.now(), 'bk cluster assignments')

        cluster01_new, freq01_new, cluster01_assignment_new, item_baskets01_new = clusters_assignment(
            baskets, item_baskets, distances_dict, nbaskets, nitems)

        if len(freq01_new[0]) == 0:
            m0 = m1_new
            m1 = None
            cluster01[0] = cluster01[1]
            cluster01[1] = None
            d0 = distances_dict[-2]
            d1 = None
            break

        if len(freq01_new[1]) == 0:
            m0 = m0_new
            m1 = None
            cluster01[0] = cluster01[0]
            cluster01[1] = None
            d0 = distances_dict[-1]
            d1 = None
            break

        cluster01 = cluster01_new
        freq01 = freq01_new
        item_baskets01 = item_baskets01_new

        d0 = distances_dict[-1]
        d1 = distances_dict[-2]

    if verbose:
        print('\t', datetime.datetime.now(), 'bk convergence')

    if d1 is not None:
        for b in cluster01[0]:
            del d1[b]

    if d0 is not None:
        if cluster01[1] is not None:
            for b in cluster01[1]:
                del d0[b]

    bisecting_kmeans_res = {
        'm0': m0,
        'm1': m1,
        'c0': cluster01[0],
        'c1': cluster01[1],
        'd0': d0,
        'd1': None if d1 is None else d1,
        'i': count
    }

    return bisecting_kmeans_res


def get_random_m0m1(baskets_id):
    m0_index = random.choice(list(baskets_id))
    m1_index = random.choice(list(baskets_id))
    while m0_index == m1_index:
        m1_index = random.choice(list(baskets_id))
    return m0_index, m1_index


def random_init(baskets, neighbors, use_neighbors, num=100):

    max_dist = -float('infinity')
    best_res = (None, None)

    for i in range(0, num):
        m0_index, m1_index = get_random_m0m1(baskets.keys())

        m0 = baskets[m0_index]
        m1 = baskets[m1_index]

        nbr_m0 = None
        nbr_m1 = None
        if use_neighbors:
            nbr_m0 = neighbors[m0_index]
            nbr_m1 = neighbors[m1_index]

        dist = calculate_distance(m0, m1, nbr_m0, nbr_m1, use_neighbors)

        if dist > max_dist:
            max_dist = dist
            best_res = (m0_index, m1_index)

    return best_res


def no_nbr_init(max_no_nbr, baskets, neighbors, use_neighbors):
    distances_tuple = dict()
    for m0_index, m1_list in max_no_nbr:
        m0 = baskets[m0_index]
        nbr_m0 = neighbors[m0_index]
        for index in range(0, len(m1_list)):
            if m1_list[index]:
                m1_index = index
                m1 = baskets[m1_index]
                nbr_m1 = neighbors[m1_index]
                distances_tuple[(m0_index, m1_index)] = calculate_distance(
                        m0, m1, nbr_m0, nbr_m1, use_neighbors)
    max_dist_couple = max(distances_tuple.items(), key=lambda x: x[1])

    return max_dist_couple[0]


def calculate_initial_centroids_2_furthest(baskets, verbose=False):
    m_rnd_index = random.choice(baskets.keys())
    m_rnd = baskets[m_rnd_index]

    m0_best = None
    max_dist = -float('infinity')
    for bid in baskets:

        m0 = baskets[bid]
        dist = calculate_distance(m0, m_rnd, None, None, False)

        if dist > max_dist:
            max_dist = dist
            m0_best = bid

    m0 = baskets[m0_best]

    m1_best = None
    max_dist = -float('infinity')
    for bid in baskets:

        m1 = baskets[bid]
        dist = calculate_distance(m0, m1, None, None, False)

        if dist > max_dist:
            max_dist = dist
            m1_best = bid

    return m0_best, m1_best


def calculate_initial_centroids_first_run(
        baskets, neighbors, no_neighbors, item_baskets, use_neighbors, num_rnd_init, verbose=False):

    max_no_nbr = list()
    max_no_nbr_index = 0

    if use_neighbors:
        for b in baskets:
            no_nbr_list = no_neighbors[b]
            no_nbr_tuple = (b, no_nbr_list)
            no_nbr = no_nbr_list.count()
            if no_nbr > max_no_nbr_index:
                max_no_nbr_index = no_nbr
                max_no_nbr = list()
                max_no_nbr.append(no_nbr_tuple)
            elif no_nbr == max_no_nbr_index:
                max_no_nbr.append(no_nbr_tuple)

    if max_no_nbr_index == 0 or len(max_no_nbr) == len(baskets):
        # print('rnd')
        if verbose:
            print(datetime.datetime.now(), 'init random')

        num_rnd_init = min(num_rnd_init, int(binom(len(baskets), 2)))
        m0_index, m1_index = random_init(baskets, neighbors, use_neighbors, num=num_rnd_init)
    else:
        # print('no nbr')
        if verbose:
            print(datetime.datetime.now(), 'init no neighbors')

        m0_index, m1_index = no_nbr_init(max_no_nbr, baskets, neighbors, use_neighbors)

    return m0_index, m1_index


def calculate_initial_centroids_max_dist(baskets, dist_from_rep, neighbors, use_neighbors, num_init_dist):

    max_dist_dict = dict()
    sorted_dist_from_rep = sorted(dist_from_rep, key=dist_from_rep.get, reverse=True)

    num_init_dist = min(num_init_dist, len(sorted_dist_from_rep))

    for m0_index in sorted_dist_from_rep[0:num_init_dist]:
        nbr_m0 = None
        for m1_index in sorted_dist_from_rep:
            if m1_index == m0_index:
                continue
            nbr_m1 = None
            if use_neighbors:
                nbr_m0 = neighbors[m0_index]
                nbr_m1 = neighbors[m1_index]
            dist = calculate_distance(baskets[m0_index], baskets[m1_index], nbr_m0, nbr_m1, use_neighbors)
            max_dist_dict[(m0_index, m1_index)] = dist

    max_dist_couple = max(max_dist_dict.items(), key=lambda x: x[1])

    return max_dist_couple[0][0], max_dist_couple[0][1]


def filter_items(baskets, item_baskets, nitems, items_in_every_basket=None):

    if items_in_every_basket is None:
        items_in_every_basket = nitems * bitarray('0')

    items_list = item_baskets.keys()
    for item in list(items_list):
        if item_baskets[item].count() >= len(baskets) * 1.0:
            items_in_every_basket[item] = 1
            del item_baskets[item]

    if not items_in_every_basket.any():
        return baskets, items_in_every_basket, item_baskets, False

    filtered_baskets = dict()
    is_homogeneous = True
    for b in baskets:
        filtered_baskets[b] = (baskets[b] ^ items_in_every_basket) & baskets[b]
        is_homogeneous = is_homogeneous and not filtered_baskets[b].any()

    return filtered_baskets, items_in_every_basket, item_baskets, is_homogeneous


def build_cluster_split(bisective_kmeans_res):
    clusters_split = [bisective_kmeans_res['c0']]
    centroids_split = [bisective_kmeans_res['m0']]
    distances_split = [bisective_kmeans_res['d0']]

    if bisective_kmeans_res['c1'] is not None:
        clusters_split.append(bisective_kmeans_res['c1'])
        centroids_split.append(bisective_kmeans_res['m1'])
        distances_split.append(bisective_kmeans_res['d1'])

    return clusters_split, centroids_split, distances_split


def build_cluster(cluster, centroid, items_in_every_basket):

    for b in cluster:
        cluster[b] |= items_in_every_basket

    centroid |= items_in_every_basket

    return {
        'cluster': cluster,
        'centroid': centroid
    }


class TXmeans:

    def __init__(self):

        self.clustering = list()
        self.medioids = list()
        self.iter_count = 0
        self.iter_bk_count = list()
        self.stack = list()

    def fit(self, baskets, nbaskets, nitems,
            num_rnd_init=100,
            min_cluster_size=3,
            min_item_freq=2,
            num_init_dist=0,
            use_neighbors=False,
            force_first_split=False,
            merge_clusters=False,
            random_sample=float('infinity'),
            verbose=False,
            random_state=None):

        self.nbaskets = nbaskets
        self.nitems = nitems

        self.num_rnd_init = num_rnd_init
        self.min_cluster_size = min_cluster_size
        self.min_item_freq = min_item_freq
        self.num_init_dist = num_init_dist
        self.use_neighbors = use_neighbors
        self.force_first_split = force_first_split
        self.merge_clusters = merge_clusters
        self.random_sample = random_sample
        self.verbose = verbose

        self.clustering = list()
        self.medioids = list()
        self.iter_count = 0
        self.iter_bk_count = list()
        self.stack = list()

        # For reproducibility
        random.seed(random_state)

        if self.verbose:
            print(datetime.datetime.now(), 'initialization')
        self._first_iter(baskets)

        while len(self.stack) > 0:

            if self.verbose:
                print(datetime.datetime.now(), 'iteration', self.iter_count)

            self.iter_count += 1
            self._iteration()

        if self.merge_clusters:
            self._merge_clusters()

        if self.nbaskets > self.random_sample:
            self._assign_baskets_to_centroids(baskets)

        self._find_medioids()

        return self


    def _first_iter(self, baskets):

        if self.nbaskets > self.random_sample:
            sample_baskets_keys = set(random.sample(baskets.keys(), self.random_sample))
            sample_baskets = dict()
            for b in sample_baskets_keys:
                sample_baskets[b] = baskets[b]
            baskets = sample_baskets
            if self.verbose:
                print(datetime.datetime.now(), 'calculate item baskets')
            item_baskets = calculate_item_baskets(baskets, self.nbaskets)
        else:
            if self.verbose:
                print(datetime.datetime.now(), 'calculate item baskets')
            item_baskets = calculate_item_baskets(baskets, self.nbaskets)

        if self.verbose:
            print(datetime.datetime.now(), 'tot baskets', self.nbaskets, 'tot item', self.nitems)

        if self.verbose:
            print(datetime.datetime.now(), 'filter items')

        baskets, items_in_every_basket, item_baskets, is_homogeneous = \
            filter_items(baskets, item_baskets, self.nitems)

        if is_homogeneous:

            if self.verbose:
                print(datetime.datetime.now(), 'the cluster is homogeneous')

            cluster_out = build_cluster(baskets, self.nitems * bitarray('0'), items_in_every_basket)
            self.clustering.append(cluster_out)
            return self

        neighbors = None
        no_neighbors = None

        if self.use_neighbors:

            if self.verbose:
                print(datetime.datetime.now(), 'calculate neighbors')

            neighbors, no_neighbors = calculate_neighbors(baskets, item_baskets, self.nbaskets, self.nitems)

        if self.verbose:
            print(datetime.datetime.now(), 'calculate frequencies')

        freq = calculate_frequencies(item_baskets)

        if self.verbose:
            print(datetime.datetime.now(), 'calculate rep')

        baskets_rep, baskets_rep_nbr = rep(
                baskets, item_baskets, freq, self.nbaskets, self.nitems, self.min_item_freq)

        if self.force_first_split:

            dist_from_rep = None
            if self.num_init_dist > 0:
                dist_from_rep = dict()
                for b in baskets:
                    nbr_b = None
                    if self.use_neighbors:
                        nbr_b = neighbors[b]

                    dist_from_rep[b] = calculate_distance(
                            baskets[b], baskets_rep, nbr_b, baskets_rep_nbr, self.use_neighbors)

            if self.verbose:
                print(datetime.datetime.now(), 'run bisecting')

            clusters_split, centroids_split, distances_split = self._bisect(
                baskets, neighbors, no_neighbors, item_baskets, dist_from_rep)

            if len(clusters_split) < 2:

                if self.verbose:
                    print(datetime.datetime.now(), 'the cluster cannot be split anymore')

                cluster_out = build_cluster(baskets, baskets_rep, items_in_every_basket)
                self.clustering.append(cluster_out)
                return

            for cluster, centroid, distances in zip(clusters_split, centroids_split, distances_split):

                if len(cluster) >= self.min_cluster_size * 2:
                    self.stack.append({
                        'b': cluster,
                        'r': centroid,
                        'd': distances,
                        'i': items_in_every_basket
                    })
                else:
                    cluster_out = build_cluster(cluster, centroid, items_in_every_basket)
                    self.clustering.append(cluster_out)

        else:

            if self.verbose:
                print(datetime.datetime.now(), 'calculate distances from rep')

            distances = dict()
            for b in baskets:
                nbr_basket = None
                if self.use_neighbors:
                    nbr_basket = neighbors[b]
                distances[b] = calculate_distance(baskets_rep, baskets[b], baskets_rep_nbr, nbr_basket,
                                                  self.use_neighbors)

            if self.verbose:
                print(datetime.datetime.now(), 'initialize stack')

            self.stack.append({
                'b': baskets,
                'r': baskets_rep,
                'd': distances,
                'i': items_in_every_basket,
                'f': (item_baskets, neighbors, no_neighbors)
            })

    def _bisect(self, baskets, neighbors, no_neighbors, item_baskets, dist_from_rep=None):

        if self.verbose:
            print(datetime.datetime.now(), 'calculate initial centroids')

        if self.num_init_dist > 0 and dist_from_rep is not None:
            m0_index, m1_index = calculate_initial_centroids_max_dist(
                    baskets, dist_from_rep, neighbors, self.use_neighbors, self.num_init_dist)
        else:
            m0_index, m1_index = calculate_initial_centroids_first_run(
                    baskets, neighbors, no_neighbors, item_baskets,
                    self.use_neighbors, self.num_rnd_init, self.verbose)

        # m0_index, m1_index = calculate_initial_centroids_2_furthest(baskets, self.verbose)

        m0 = baskets[m0_index]
        m1 = baskets[m1_index]

        nbr_m0 = None
        nbr_m1 = None
        if self.use_neighbors:
            nbr_m0 = neighbors[m0_index]
            nbr_m1 = neighbors[m1_index]

        if self.verbose:
            print(datetime.datetime.now(), 'calculate distances from initial centroids')

        distances_dict_init = calculate_distances(
                baskets, [neighbors, neighbors], m0, m1, nbr_m0, nbr_m1, self.use_neighbors)

        if self.verbose:
            print(datetime.datetime.now(), 'run bisecting kmeans')

        bk_res = bisecting_kmeans(baskets, neighbors, item_baskets, distances_dict_init,
                                  self.nbaskets, self.nitems, m0, m1, self.use_neighbors,
                                  self.min_item_freq, self.verbose)

        self.iter_bk_count.append(bk_res['i'])

        if self.verbose:
            print(datetime.datetime.now(), 'build cluster split')

        clusters_split, centroids_split, distances_split = build_cluster_split(bk_res)

        return clusters_split, centroids_split, distances_split

    def _iteration(self):

        if self.verbose:
            print(datetime.datetime.now(), 'get element from the stack')

        elem = self.stack.pop()
        baskets = elem['b']
        baskets_rep = elem['r']
        distances = elem['d']
        items_in_every_basket = elem['i']

        if 'f' in elem:

            if self.verbose:
                print(datetime.datetime.now(), 'get precomputed item_baskets, neighbors, no_neighbors')

            item_baskets, neighbors, no_neighbors = elem['f']

        else:

            if self.verbose:
                print(datetime.datetime.now(), 'calculate item baskets')

            item_baskets = calculate_item_baskets(baskets, self.nbaskets)

            if self.verbose:
                print(datetime.datetime.now(), 'filter items')

            baskets, items_in_every_basket, item_baskets, is_homogeneous = filter_items(
                    baskets, item_baskets, self.nitems, items_in_every_basket)

            if is_homogeneous:

                if self.verbose:
                    print(datetime.datetime.now(), 'the cluster is homogeneous')

                cluster_out = build_cluster(baskets, baskets_rep, items_in_every_basket)
                self.clustering.append(cluster_out)
                return

            neighbors = None
            no_neighbors = None

            if self.use_neighbors:

                if self.verbose:
                    print(datetime.datetime.now(), 'calculate neighbors')

                neighbors, no_neighbors = calculate_neighbors(
                        baskets, item_baskets, self.nbaskets, self.nitems)

        if self.verbose:
            print(datetime.datetime.now(), 'count current baksets and items')

        ndim = 0
        for item in item_baskets:
            if item_baskets[item].any():
                ndim += 1

        npoints = len(baskets)

        if self.verbose:
            print(datetime.datetime.now(), 'baskets', npoints, 'items', ndim)

        if self.verbose:
            print(datetime.datetime.now(), 'calculate bic merge')

        bic_merge = bic([baskets], npoints, ndim, list(distances.values()))

        if self.verbose:
            print(datetime.datetime.now(), 'run bisecting')

        clusters_split, centroids_split, distances_split = self._bisect(
                baskets, neighbors, no_neighbors, item_baskets, distances)

        if len(clusters_split) < 2:

            if self.verbose:
                print(datetime.datetime.now(), 'the cluster cannot be split anymore')

            cluster_out = build_cluster(baskets, baskets_rep, items_in_every_basket)
            self.clustering.append(cluster_out)
            return

        if self.verbose:
            print(datetime.datetime.now(), 'calculate bic split')

        bic_split = bic(clusters_split, npoints, ndim, [list(distances_split[0].values()),
                                                        list(distances_split[1].values())])

        if self.verbose:
            print(datetime.datetime.now(), 'bic merge', bic_merge, 'bic split', bic_split)

        # print('merge', bic_merge, 'split', bic_split, abs(bic_split-bic_merge), \)
        #     len(clusters_split[0]), len(clusters_split[1]), bic_split > bic_merge \
        #         and len(clusters_split[0]) >= self.min_cluster_size and len(clusters_split[1]) >=
        # self.min_cluster_size

        if bic_split > bic_merge \
                and len(clusters_split[0]) >= self.min_cluster_size \
                and len(clusters_split[1]) >= self.min_cluster_size:

            if self.verbose:
                print(datetime.datetime.now(), 'another split is required', \
                    len(clusters_split[0]), len(clusters_split[1]))

            for cluster, centroid, distances in zip(clusters_split, centroids_split, distances_split):

                if len(cluster) >= self.min_cluster_size * 2:
                    self.stack.append({
                        'b': cluster,
                        'r': centroid,
                        'd': distances,
                        'i': items_in_every_basket
                    })
                else:
                    cluster_out = build_cluster(cluster, centroid, items_in_every_basket)
                    self.clustering.append(cluster_out)
        else:

            if self.verbose:
                print(datetime.datetime.now(), 'the cluster is considered enough homogeneous')

            cluster_out = build_cluster(baskets, baskets_rep, items_in_every_basket)
            self.clustering.append(cluster_out)

        if self.verbose:
            print('')

    def _merge_clusters(self):

        if self.verbose:
            print(datetime.datetime.now(), 'merge clusters')

        clusters_to_be_merged = defaultdict(list)
        for i, cluster in enumerate(self.clustering):
            centroid0 = cluster['centroid']
            for j, cluster_iter in enumerate(self.clustering):
                centroid1 = cluster_iter['centroid']
                if i != j and centroid0 == centroid1:
                    clusters_to_be_merged[i].append(j)

        clustering_new = list()
        for i, cluster in enumerate(self.clustering):
            if i not in clusters_to_be_merged:
                clustering_new.append(cluster)

        clusters_already_merged = set()
        for i in clusters_to_be_merged:
            if i not in clusters_already_merged:
                centroid0 = self.clustering[i]['centroid']
                cluster0 = self.clustering[i]['cluster']
                for j in clusters_to_be_merged[i]:
                    if j not in clusters_already_merged:
                        cluster1 = self.clustering[j]['cluster']
                        for b, basket in cluster1.iteritems():
                            cluster0[b] = basket
                        clusters_already_merged.add(j)
                clusters_already_merged.add(i)
                clustering_new.append({
                    'cluster': cluster0,
                    'centroid': centroid0
                })
        self.clustering = clustering_new

    def _assign_baskets_to_centroids(self, baskets):

        if self.verbose:
            print(datetime.datetime.now(), 'assign baskets to centroids')

        for b in baskets:
            min_dist = float('infinity')
            best_cluster = None
            for i, cluster in enumerate(self.clustering):
                centroid = cluster['centroid']
                dist = calculate_distance(baskets[b], centroid, None, None, False)
                if dist < min_dist:
                    best_cluster = i
                    min_dist = dist
            self.clustering[best_cluster]['cluster'][b] = baskets[b]

    def _find_medioids(self):
        res = self.clustering
        for label, cluster in enumerate(res):
            d = np.infty
            idx_med = 0
            for bid, bitarr in cluster['cluster'].items():
                d_tmp = jaccard_bitarray(bitarr, cluster['centroid'])
                if d_tmp < d:
                    d = d_tmp
                    idx_med = bid
            self.medioids.append(idx_med)
