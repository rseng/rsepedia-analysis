[latest release]:(https://github.com/jGetMove/jGetMove/releases/latest)

jGetMove v2.0.2 - Xénophon
============================

GetMove algorithm in java

Requirements
------------

Works on Linux, MacOS and Windows with Java 1.8 and higher

Sources
-------

To compile the project, simply use `make`. It will create `jGetMove.jar`.

If you do not have `make` installed or want build it manually use :

```bash
$ mkdir -p out/
$ javac -extdirs lib/ -sourcepath src/ src/fr/jgetmove/jgetmove/Main.java -d out/
$ jar cvfm jGetMove.jar Manifest.mf -C out/ . lib/
```

or download it here: [Latest Release]

Usage 
-----

To execute the program (e.g. to see the usage): 
```bash
$ java -jar jGetMove.jar --help
```

Input files
-----------

### Transaction-Clusters

The transaction cluster association is defined in this file (named `index.dat` by conventions).
the `Transaction` is the line number (begginning by `0`) and it iterates through the `Cluster`s written on it's line.

#### Example

> In this example there are two transactions (`0` and `1`) which respectively iterate through `[0,1,2]` and `[2,3,4]`.
> Following this logic there are 5 clusters (from `0` to `4`).
>
> **Note :** the javadoc refers to this structure with : `transactionId [clusterId ...]`

```
0 1 2
2 3 4
```
> A working example is present in `assets/example.dat`

### Cluster-Time

The relation between the cluster and the Time are defined here (named `time_index.dat` by conventions).

Each line of the file is a link between the `Time` and the `Cluster` by their id.

#### Example

> The following file is the cluster-time link of cluster-transaction file
> each line has a `Time` (beginning at `1`) and the clusters.
> The first and the second `Time` have two clusters.
> Please note that the cluster **NEED TO BE** present in `index.dat`, an exception will be returned otherwise. 

```
1 0
1 1
2 2
2 3
3 4
```
> A working example can be found in `assets/example_time_index.dat`

Reference
---------

__JGetMove: Mining Multiple Movement Patterns__

*Fati Chen, Nhathai Phan, Pascal Poncelet, Maguelonne Teisseire*

> Abstract :
> Recent improvements in positioning technology have led to a much wider availability of massive moving object data. A crucial task is to find the moving ob jects that travel together. Usually, these ob ject sets are called object movement patterns. Due to the emergence of many different kinds of object movement patterns in recent years, different approaches have been proposed to extract them. However, each approach only focuses on mining a specific kind of patterns. In addition to being a painstaking task due to the large number of algorithms used to mine and manage patterns, it is also time consuming. Moreover, we have to execute these algorithms again whenever new data are added to the existing database. To address these issues, we first redefine movement patterns in the itemset context. Secondly, we propose a unifying approach, named GeTMove, which uses a frequent closed itemset-based object movement pattern-mining algorithm to mine and manage different patterns. 
> 
> This code is a Java implementation of GetMove called JGetMove.

The description of the approach is available at :

Nhathai Phan, Pascal Poncelet, Maguelonne Teisseire. *All in One: Mining Multiple Movement Patterns. International Journal of Information Technology and Decision Making*, World Scientific Publishing, 2016, 15 (5), pp.1115-1156.
>
> Keywords: Object movement pattern; frequent closed itemset; unifying approach; trajectories.


License
-------

[![Licence Creative Commons](https://i.creativecommons.org/l/by-nc-sa/3.0/fr/88x31.png)](http://creativecommons.org/licenses/by-nc-sa/3.0/fr/)  
Ce(tte) œuvre est mise à disposition selon les termes de la [Licence Creative Commons Attribution - Pas d’Utilisation Commerciale - Partage dans les Mêmes Conditions 3.0 France](http://creativecommons.org/licenses/by-nc-sa/3.0/fr/).

For more informations see [LICENSE.md](LICENSE.md)

