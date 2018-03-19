													################################
													# Lignes de commande pour Git: #
													################################
													
													Ceci se fait en local sur le PC
													
													
Se mettre dans le bon dossier (celui du repository ou celui dans lequel on veut cloner le repository) dans l'invite de commande.

git clone clone https://gitlab.icm-institute.org/fanny.grosselin/CARE-rCortex.git  --> pour cloner en local le repository 
OU
git pull origin nom-de-la-branche   --> pour récupérer le repo de Git depuis nom-de-la-branche (quand on possède déjà le repo en local mais qu'on veut les dernières modifications)


git branch    --> pour connaître les branches existantes
git checkout nom-de-la-branche   --> pour changer de branche

POUR METTRE LES MODIFICATIONS SUR GIT:
-------------------------------------
1) git status     --> dit les fichiers qui diffèrent entre le dossier local (sur l'ordinateur) et ce qui est sur le serveur

2) git add .    --> pour ajouter avant le commit les fichier qu'on veux commiter
	OU
   git add nom-du-fichier   --> pour ajouter des fichiers parmis ceux du status qu'on veut commiter

3) git commit -m "Description du commit"  --> pour commiter

4) git push -u origin nom-de-la-branche --> pour push sur Git sur nom-de-la-branche



_____________________________________________________________________________________________________________________________________
Pour ne pas tracker certains fichiers (qu'on ne veut pas mettre sur Git comme par exemple des fichiers test):
Créer en local un fichier texte avec pour nom ".gitignore" (ne pas oublier le .)
Dans ce fichier, il suffit de mettre sur chaque ligne l'adresse complète du fichier ou du dossier qu'on veut que Git ignore.
Si on fait des changements dans ces fichiers ou dossiers, lorsqu'on fait git status, il ne seront pas affiché.

______________________________________________________________________________________________________________________________________
Lien utile pour trouver les lignes de commande pour git : https://gist.github.com/aquelito/8596717