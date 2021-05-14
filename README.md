# TX-Means

TX-Means is a parameter-free clustering algorithm able to efficiently partitioning transactional data in a completely automatic way.
TX-Means is designed for the case where clustering must be applied on a massive number of different datasets, for instance when a large set of users need to be analyzed individually and each of them has generated a long history of transactions.

In this repository we provide the source code of TX-Means, the clustering algorithm competitors and the dataset used in
> Riccardo Guidotti, Anna Monreale, Mirco Nanni, Fosca Giannotti, Dino Pedreschi *"Clustering Individual Transactional Data for Masses of Users"*, KDD 2017, 2017, Halifax, NS, Canada

Please cite the paper above if you use our code or dataets.

Requirements:
- python >= 3 
- numpy >= 1.10.1
- pandas >= 0.18.1
- scipy >= 0.17.1
- bitarray >= 0.8.1
- Java >= 8.1
