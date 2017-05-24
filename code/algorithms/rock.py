import copy
import random

from heapq import *

from util import *

__author__ = 'Riccardo Guidotti'


def calculate_links(baskets, neighbors):
    links_dict = dict()

    for b1 in baskets:
        links_dict[b1] = dict()
        for b2 in baskets:
            links = 1.0 * (neighbors[b1] & neighbors[b2]).count()
            links_dict[b1][b2] = links

    return links_dict


def goodness(b1, b2, c1, c2, links_dict, sigma):
    num = links_dict[b1][b2]
    den = (len(c1) + len(c2))**(1 + 2*fun(sigma)) - (len(c1))**(1 + 2*fun(sigma)) - (len(c2))**(1 + 2*fun(sigma))
    return 1.0 * num / den


def fun(sigma):
    return 1.0 * (1 - sigma) / (1 + sigma)


def calculate_neighbors(baskets, nbaskets, nitems, sigma):
    neighbors = defaultdict(lambda: nbaskets * bitarray('0'))
    for b in baskets:
        for b1 in baskets:
            if b == b1:
                continue
            if jaccard_bitarray(baskets[b], baskets[b1]) >= sigma:
                neighbors[b][b1] = 1
        neighbors[b][b] = 0

    return neighbors


def sample_size(population, conf_interval, conf_level=0.95, prob=0.5):
    if conf_interval < 0.01:
        return population

    zscore = stats.norm.ppf(1-(1-conf_level)/2)
    ss = ((zscore ** 2) * prob * (1-prob)) / (conf_interval ** 2)
    new_ss = ss / (1 + (ss - 1) / population)

    return int(round(new_ss))


class Rock:

    def __init__(self):
        self.clustering = list()
        self.iter_count = 0

    def fit(self, baskets, nbaskets, nitems, k, sigma, random_sample=float('infinity')):

        self.clustering = list()
        self.iter_count = 0

        self.nbaskets = nbaskets
        self.nitems = nitems

        self.k = k
        self.sigma = sigma

        self.random_sample = random_sample

        clusters = self._run(baskets)

        if self.nbaskets > self.random_sample:
            # print 'labeling baskets'
            clusters = self._labeling_baskets(baskets, clusters)

        for cluster in clusters.values():
            res = dict()
            for c in cluster:
                res[c] = baskets[c]
            self.clustering.append({
                'cluster': res,
                'centroid': None
            })

        return self

    def _run(self, baskets):

        if self.nbaskets > self.random_sample:
            sample_baskets_keys = set(random.sample(baskets.keys(), self.random_sample))
            sample_baskets = dict()
            for b in sample_baskets_keys:
                sample_baskets[b] = baskets[b]
            baskets = sample_baskets

        clusters = dict()
        for b in baskets:
            clusters[b] = [b]

        prog_cluster_id = self.nbaskets

        # print 'calculating neighbors'
        neighbors = calculate_neighbors(baskets, self.nbaskets, self.nitems, self.sigma)

        # print 'calculating links'
        links_dict = calculate_links(baskets, neighbors)

        q = dict()

        for b1 in links_dict:
            q[b1] = list()
            for b2 in links_dict[b1]:
                if b2 == b1:
                    continue

                val = -goodness(b1, b2, clusters[b1], clusters[b2], links_dict, self.sigma)
                heappush(q[b1], (val, b2))

        max_dict = dict()

        for b in q:
            if len(q[b]) > 0:
                max_dict[b] = q[b][0][1]

        Q = list()
        for b in q:
            if len(q[b]) > 0:
                heappush(Q, (q[max_dict[b]][0][0], b))

        self.iter_count +=1

        while len(Q) > self.k:
            # print 'iter', self.iter_count, len(Q)
            self.iter_count += 1
            u = heappop(Q)[1]
            v = max_dict[u]

            for e in Q:
                if e[1] == v:
                    Q.remove(e)

            w = prog_cluster_id
            prog_cluster_id += 1

            # print 'u,v,w', u, v, w

            q[w] = list()
            clusters[w] = list()
            clusters[w].extend(copy.deepcopy(clusters[u]))
            clusters[w].extend(copy.deepcopy(clusters[v]))
            links_dict[w] = dict()

            # print u, q[u]
            # print v, q[v]

            in_qu = set()
            for xc in q[u]:
                x = xc[1]
                in_qu.add(x)

            in_qv = set()
            for xc in q[v]:
                x = xc[1]
                in_qv.add(x)

            for x in in_qu:
                # print 'x', x
                if x == u or x == v:
                    continue

                lxu = 0 if u not in links_dict[x] else links_dict[x][u]
                lxv = 0 if v not in links_dict[x] else links_dict[x][v]
                links_dict[x][w] = lxu + lxv
                links_dict[w][x] = lxu + lxv

                new_max = float('infinity')
                new_max_elem = None

                to_be_removed = set()
                for e in q[x]:
                    if e[1] == u:
                        to_be_removed.add(e)
                    elif x in in_qv and e[1] == v:
                        to_be_removed.add(e)
                    elif new_max > e[0]:
                        new_max = e[0]
                        new_max_elem = e[1]

                for tbr in to_be_removed:
                    q[x].remove(tbr)

                max_dict[x] = new_max_elem if new_max_elem is not None else max_dict[x]
                val = -goodness(x, w, clusters[x], clusters[w], links_dict, self.sigma)
                heappush(q[x], (val, w))
                heappush(q[w], (val, x))

                if q[max_dict[x]][0][0] > val:
                    max_dict[x] = w
                    for e in Q:
                        if e[1] == x:
                            Q.remove(e)

                    heappush(Q, (val, x))

            for x in in_qv:
                if x not in in_qu:
                    if x == u or x == v or x >= w:
                        continue

                    # print 'x2', x
                    lxu = 0 if u not in links_dict[x] else links_dict[x][u]
                    lxv = 0 if v not in links_dict[x] else links_dict[x][v]
                    links_dict[x][w] = lxu + lxv
                    links_dict[w][x] = lxu + lxv
                    new_max = float('infinity')
                    new_max_elem = None

                    to_be_removed = set()
                    for e in q[x]:
                        if e[1] == v:
                            to_be_removed.add(e)
                        elif new_max > e[0]:
                            new_max = e[0]
                            new_max_elem = e[1]

                    for tbr in to_be_removed:
                        q[x].remove(tbr)

                    max_dict[x] = new_max_elem if new_max_elem is not None else max_dict[x]
                    val = -goodness(x, w, clusters[x], clusters[w], links_dict, self.sigma)
                    heappush(q[x], (val, w))
                    heappush(q[w], (val, x))

                    if q[max_dict[x]][0][0] > val:
                        max_dict[x] = w
                        for e in Q:
                            if e[1] == x:
                                Q.remove(e)
                        heappush(Q, (val, x))

            max_dict[w] = q[w][0][1]
            heappush(Q, (q[max_dict[w]][0][0], w))

            del q[u]
            del q[v]
            del max_dict[u]
            del max_dict[v]
            del clusters[u]
            del clusters[v]

        return clusters

    def _labeling_baskets(self, baskets, clusters):

        clusters_sample = dict()

        for cid, cluster in clusters.iteritems():
            if len(cluster) > 1000:
                ss = sample_size(len(cluster), 0.05, conf_level=0.95, prob=0.5)
                clusters_sample[cid] = random.sample(cluster, ss)
            else:
                clusters_sample[cid] = cluster

        clusters_new = defaultdict(list)
        for bid, basket in baskets.iteritems():

            best_cid = None
            max_nbr_count = -float('infinity')
            for cid in clusters_sample:
                nbr_count = 0
                for binc in clusters_sample[cid]:
                    if jaccard_bitarray(basket, baskets[binc]) >= self.sigma:
                        nbr_count += 1
                if nbr_count > max_nbr_count:
                    max_nbr_count = nbr_count
                    best_cid = cid

            clusters_new[best_cid].append(bid)

        return clusters_new
