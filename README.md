# projet GSA avec PBxplore : Analyse de la communication allostérique à l'aide d'un alphabet

Dai Yuping 

14/09/2021




## 1.Installer les packages :

fichier env_projet_court.yml contenant les différents prérequis et pouvant être exploité par conda pour créer un environnement viable pour le programme.

conda env create -f env_projet_court.yml

Environement :

pbxplore 1.4.0

numpy  1.21.2 

pandas 1.3.3

matplotlib  3.4.3  

seaborn 0.11.2


## 2.Étape PBxplore
La première étape est utilisée PBassign pour le fichier des PB :

[git_PBxplore](https://github.com/pierrepo/PBxplore)

PBassign -p Frames/ -o test

-p Frames/ : Entrée un dossier contenant des fichiers PDB 

-o test : le nom de fichiers fasta (test.PB.fasta)

Output : fichier fasta (test.PB.fasta) qui va utiliser comme le dossier entrés pour la prochaine étape

## 3. Utilisation toto.py 

Le programme va utiliser le fichier fasta pour calculer les fréquences, faire le matrix fréquences et faire la matrix d’information mutuelle

python toto.py --fasta '/ path/test.PB.fasta' --vitesse_transitions 2

--fasta[path fille] : correspond au fichier fasta qui vient de l’étape dernière.

--vitesse_transitions [int] : Vitesse transition (valeur défaut : 1)

-h ou –help : permetra d’afficher l’aide à l’utilisation.

Output : mutual_information_.pdf :

              matrix_MI.csv


