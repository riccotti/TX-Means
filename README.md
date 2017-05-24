# TX-Means

TX-Means is a parameter-free clustering algorithm able to efficiently partitioning transactional data in a completely automatic way.
TX-Means is designed for the case where clustering must be applied on a massive number of different datasets, for instance when a large set of users need to be analyzed individually and each of them has generated a long history of transactions.

In this repository we provide the source code of TX-Means, the clustering algorithm competitors and the dataset used in
> Riccardo Guidotti, Anna Monreale, Mirco Nanni, Fosca Giannotti, Dino Pedreschi *"Clustering Individual Transactional Data for Masses of Users"*, KDD 2017, 2017, Halifax, NS, Canada

Please cite the paper above if you use our code or dataets.


- Files contained in the folder:
	README.txt
	code/
		__init__.py
		test.py
		algorithms/
			__init__.py
			atdc.py
			txmeans.py
			clope.py
			coolcat.py
			rock.py
			tkmeans.py
			util.py
		generators/
			__init__.py
			ibm_generator
			atdc_datagenerator.jar
			datagenerator.py
			datamanager.py
		validation/
			__init__.py
			validation_measures.py
			calculate_aggregate_statistics.py		
	dataset/
		mushrooms.csv
		congress.csv
		zoo.csv
		sample_real_data.csv
	

- Requirements:

python >= 2.7 (not python 3)
numpy >= 1.10.1
pandas >= 0.18.1
scipy >= 0.17.1
bitarray >= 0.8.1
Java >= 8.1
