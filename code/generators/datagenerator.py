import os
import numpy as np
import subprocess
import shutil
from datamanager import *

__author__ = 'Riccardo Guidotti'

def generate_syntetic_data1(nclus, ntrans, nitems, avg_tlen,
                            min_tlen=-float('infinity'), max_tlen=float('infinity'), pattern_cor_level=1.0):
    """
    nclus = 4 # number of clusters
    ntrans = 16 # number of transcations
    nitems = 32 # max number of items
    avg_tlen = 6 # average transaction length
    min_tlen = 6 # minimum transaction length
    max_tlen = 6 # maximum transaction length
    pattern_cor_level = 1.0 # level of correlation inside a cluster

    if min_tlen == max_tlen == avg_tlen all transactions have the same length
    if pattern_cor_level == 1.0 every cluster contains only items assigned to it
    if pattern_cor_level == 0.0 every cluster contains only items assigned to another"""

    if 1.0*nitems/nclus < avg_tlen:
        print 'Attention!!! nitems/nclus < avg_tlen', 1.0*nitems/nclus, avg_tlen
        return None

    # init transactions
    # divido le transazioni in nclus
    trans = range(0, ntrans)
    trans_clusters = [trans[x:x+ntrans/nclus] for x in range(0, ntrans, ntrans/nclus)]
    if len(trans_clusters) > nclus:
        for tid, index in zip(trans_clusters[nclus], range(0, len(trans_clusters[nclus]))):
            trans_clusters[index].append(tid)
        trans_clusters.pop()

    # init items
    # associo ai cluster gli item che li caratterizzano
    items = range(0, nitems)
    items_clusters = [items[x:x+nitems/nclus] for x in range(0, nitems, nitems/nclus)]
    if len(items_clusters) > nclus:
        for tid, index in zip(items_clusters[nclus], range(0, len(items_clusters[nclus]))):
            items_clusters[index].append(tid)
        items_clusters.pop()

    # per ogni gruppo
    # per ogni transazione
    # pesco la lunghezza come numero di volte estratto da una Poisson distribution con media tlen
    # pesco con prob pattern_cor_level in maniera uniforme fra gli item che caratterizzano il gruppo,
    # altrimenti tra gli altri
    baskets = list()
    cluster_label = 0
    for cluster, item_cluster in zip(trans_clusters, items_clusters):
        p_in_cluster = pattern_cor_level / len(item_cluster)
        p_out_cluster = (1.0-pattern_cor_level) / (nitems - len(item_cluster))

        p_items = [p_out_cluster] * nitems
        for tid in item_cluster:
            p_items[tid] = p_in_cluster

        for tid in cluster:
            tlen = np.max([np.min([np.random.poisson(avg_tlen), max_tlen]), min_tlen])
            basket = np.random.choice(items, tlen, replace=False, p=p_items).tolist()
            baskets.append((list(basket), cluster_label))
        cluster_label += 1

    return baskets


def generate_syntetic_data2(nclus, nsub, nnosub, ntrans, nitems, avg_tlen,
                            min_tlen=-float('infinity'), max_tlen=float('infinity'),
                            pattern_cor_level=1.0, perc_common_items_in_subcluster=0.6):
    """
    nclus = 4 # number of clusters
    nsub = 2 # number of subcluster for each cluster
    nnosub = 1 # number of cluster to not be split
    ntrans = 16 # number of transcations
    nitems = 32 # max number of items
    avg_tlen = 4 # average transaction length
    min_tlen = 4 # minimum transaction length
    max_tlen = 4 # maximum transaction length
    pattern_cor_level = 1.0 # level of correlation inside a cluster
    perc_common_items_in_subcluster = 0.6 # percentage of shared items in subclusters

    if min_tlen == max_tlen == avg_tlen all transactions have the same length
    if pattern_cor_level == 1.0 every cluster contains only items assigned to it
    if pattern_cor_level == 0.0 every cluster contains only items assigned to another
    if nclus == nnosub and nsub == 0 same effect of generate_syntetic_data1"""

    if 1.0*nitems/((nclus-nnosub)*nsub + nnosub) < avg_tlen:
        print 'Attention!!! nitems/nclus < avg_tlen', 1.0*nitems/((nclus-nnosub)*nsub + nnosub), avg_tlen
        return None

    # init transactions
    # divido le transazioni in nclus
    trans = range(0, ntrans)
    trans_clusters = [trans[x:x+ntrans/nclus] for x in range(0, ntrans, ntrans/nclus)]
    if len(trans_clusters) > nclus:
        for tid, index in zip(trans_clusters[nclus], range(0, len(trans_clusters[nclus]))):
            trans_clusters[index].append(tid)
        trans_clusters.pop()

    # divido i nclus - nnosub in nsub subclusters
    trans_sub_clusters = dict()
    cluster_label = 0
    for cluster in trans_clusters:
        clen = len(cluster)
        if cluster_label < nclus-nnosub:
            trans_sub_clusters[cluster_label] = [cluster[x:x+clen/nsub] for x in range(0, clen, clen/nsub)]
        else:
            trans_sub_clusters[cluster_label] = [cluster]
        cluster_label += 1

    # init items
    # associo ai cluster gli item che li caratterizzano
    items = range(0, nitems)
    items_clusters = [items[x:x+nitems/nclus] for x in range(0, nitems, nitems/nclus)]
    if len(items_clusters) > nclus:
        for tid, index in zip(items_clusters[nclus], range(0, len(items_clusters[nclus]))):
            items_clusters[index].append(tid)
        items_clusters.pop()

    # spartisco gli items nei subclusters
    cluster_label = 0
    items_sub_clusters = dict()
    for items_cluster in items_clusters:
        sub_clus_len = int(np.round(len(items_cluster) * perc_common_items_in_subcluster))
        if cluster_label < nclus-nnosub:
            items_sub_clusters[cluster_label] = [items_cluster[:sub_clus_len], items_cluster[-sub_clus_len:]]
        else:
            items_sub_clusters[cluster_label] = [items_cluster]
        cluster_label += 1

    # per ogni gruppo
    # per ogni transazione
    # pesco la lunghezza come numero di volte estratto da una Poisson distribution con media tlen
    # pesco con prob pattern_cor_level in maniera uniforme fra gli item che caratterizzano il gruppo,
    # altrimenti tra gli altri
    baskets = list()
    cluster_label = 0
    for clid in range(0, nclus):
        for sub_cluster, items_sub_cluster in zip(trans_sub_clusters[clid], items_sub_clusters[clid]):
            p_in_cluster = pattern_cor_level / len(items_sub_cluster)
            p_out_cluster = (1.0-pattern_cor_level) / (nitems - len(items_sub_cluster))

            p_items = [p_out_cluster] * nitems
            for tid in items_sub_cluster:
                p_items[tid] = p_in_cluster

            for tid in sub_cluster:
                tlen = np.max([np.min([np.random.poisson(avg_tlen), max_tlen]), min_tlen])
                tsublen = tlen
                if clid >= nclus-nnosub:
                    tsublen = tlen * 2
                basket = np.random.choice(items, tsublen, replace=False, p=p_items).tolist()
                baskets.append((list(basket), cluster_label))
            cluster_label += 1

    return baskets


def generate_syntetic_data3(nitems=1000, ntrans=1000, tlen=2, nclasses=2, out=60, l=50, filename='filename'):
    """
    nitems = 10 # average number of items (default 1000)
    ntrans = 10 # average number of transactions (default 1000)
    tlen = 5 # average size of transactions (default 10)
    nclasses = 2 # average number of classes (default 2)
    out = 0 # percentage of outliers (default 60)
    l = 0 # percentage of overlap (default 50)"""

    # java -Xmx1024M -jar dg_manco2.jar -10 -nitems 10 -ntrans 10 -tlen 5 -nclasses 2 -out 0 -l 0 -format plain

    path = './'

    subprocess.call(['java', '-Xmx1024m', '-jar', path + 'atdc_datagenerator.jar',
                     '-nitems', '%s' % nitems,
                     '-ntrans', '%s' % ntrans,
                     '-tlen', '%s' % tlen,
                     '-nclasses', '%s' % nclasses,
                     '-out', '%s' % out,
                     '-l', '%s' % l,
                     '-format', 'plain',
                     '-file', filename])


def generate_syntetic_data4(ntrans=1000, tlen=10, nitems=100):
    """
    -ntrans number_of_transactions_in_000s (default: 1000)
    -tlen avg_items_per_transaction (default: 10)
    -nitems number_of_different_items_in_000s) (default: 100)
    """

    # ./ibm_datagenerator lit -ntrans 0.01 -tlen 4 -nitems 0.1

    out = subprocess.check_output(['./ibm_datagenerator', 'lit',
                     '-ntrans', '%s' % ntrans,
                     '-tlen', '%s' % tlen,
                     '-nitems', '%s' % nitems])

    baskets = list()
    basket = list()
    basket_id = None
    rows = out.rstrip().split('\n')
    for row in rows:
        fields = row.rstrip().split(' ')
        if basket_id is None:
            basket_id = fields[0]
        new_basket_id = fields[0]
        if new_basket_id != basket_id:
            baskets.append((basket, -1))
            basket_id = new_basket_id
            basket = list()
        basket.append(int(fields[1]))

    os.remove('data.ntrans_%s.tlen_%s.nitems_%s' % (ntrans, tlen, nitems))
    os.remove('pat.ntrans_%s.tlen_%s.nitems_%s' % (ntrans, tlen, nitems))

    return baskets
