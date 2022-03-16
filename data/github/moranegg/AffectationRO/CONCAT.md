projet en Recherche Opérationelle
=================================

[![DOI](https://zenodo.org/badge/35734299.svg)](https://zenodo.org/badge/latestdoi/35734299)

Problème d'affectation généralisée
----------------------------------
un système distribué comme un ensemble de processeurs
pouvant exécuter des tâches (ou processus) en parallèle. On considère donc un ensemble de m
processeurs, chacun muni d’une certaine quantité de mémoire vive (RAM), qu’il peut utiliser
pour charger et exécuter des tâches, et un ensemble de n tâches à exécuter, chacune nécessitant
une certaine quantité de RAM pour être chargée et exécutée. Cette quantité peut en fait varier
en fonction de la nature du processeur sur lequel la tâche est exécutée, et dépend donc du choix
de ce processeur.
Enfin, à chaque couple (processeur, tâche), on associe un coût à payer pour exécuter cette
tâche sur ce processeur, et à chaque couple de tâches on associe un coût de communication (coût
à payer si ces deux tâches sont exécutées sur des processeurs différents).

Paramètres:

 * n : nombre de tâches
 * m : nombre de processeurs
 * b[] : quantité de RAM disponible sur chaque processeur
 * c[][] :coût d'affectation d'une tâche à un processeur
 * a[][] :RAM nécessaire pour une tâches i sur processeur j
 * nprim : nombre minimale de tâches à exécuter
 
Variables:

 * x[][] : affectation d'une tâche i à un processeur j
 * y[] : dans le cas d'une exécution d'un minimum de tâches n'
 
Contraintes:

 * quantité de Ram sur processeur j est limité à b[j]-> la somme de a[i,j]<=b[j] pour un processeur j
 * Une tâche s'exécute seulement une seule fois -> la somme de x[i,j] = 1 pour une tâche i
 
Objectif:
 
 * trouver le Min des coût d'exécution


 

Instructions (dans un système Windows)
------------

Assurer vous d'avoir GLPK installé:

 * GLPK pour windows : <a href="http://winglpk.sourceforge.net/">télécharger ici</a>
 * Configurer le chemin pour exécuter glpsol.exe 

Pour utiliser les fichiers :
* option A:
	 * ouvrez cmd et placer vous dans le dossier où se trouve le fichier
	 * exécutez et recupérez la solution: glpsol -m projet.mod -d nomFichierData.DAT > fictest.txt
* option B:
	* lancez en cmd projet.bat projet.mod nomFichierData.DAT fictest.txt
	
* Puis, pour une solution entière
 * placez fictest.txt dans le repertoire du programme java
 * lancez programme java pour afficher une solution entière avec fictest.txt


Fichiers important:

| Fichier | Description |
| ---- | ----------- |
| projet.bat | exécution de la solution sur windows |
| projet.mod | modélisation du problème avec obligation d'exécution de toutes les tâches |
| projetnprim.mod | exécuter au minimum n' tâches |


classes Importante:

| Class | Description |
| ----- | ----------- |
| Client | main-> affichage de la solution à partir d'un fichier texte|
| FichierSol | Lecture d'un fichier texte et transformation en Solution |
| Solution | Objet solution - gère la transformation en solution entière |



Histoire des versions
---------------

 * v1.0 première étape du projet


