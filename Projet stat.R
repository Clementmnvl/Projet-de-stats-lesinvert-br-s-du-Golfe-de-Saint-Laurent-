#Importation des biblioth�ques utilis�es
library(FactoMineR)
library(car)
library(emmeans)
library(corrplot)
library(ppcor)
library(knitr)

#Imporation de nos donn�es � �tudier 
donnees <- read.table("StLaurent.csv",                    # Chemin vers le fichier
                       sep = ";",                         # S�parateur de champs
                       header = TRUE,                     # La 1ere ligne est le nom de colonnes
                       colClasses = c(year ="factor", strate = "factor", sediment ="factor", starfish ="factor", urchin="factor"), # On pr�cise le type pour les colonnes ambigues
                       skip = 0,                          # On ne saute aucune lignes
                       dec = ",")                         # Le s�pareteur d�cimal est la ,

summary(donnees) #On obtient un resume du jeu de donnees
str(donnees)

#Representation de la variable invertebrate en fonction des autres variables:
plot(donnees$depth,donnees$invertebrate,xlab='profondeur',ylab='biomasse d invert�br�s') 
plot(donnees$totconsum,donnees$invertebrate,xlab='pr�dateur d invert�br�s',ylab='biomasse d invert�br�s') 
plot(donnees$temperature,donnees$invertebrate,xlab='temperature',ylab='biomasse d invert�br�s')
plot(donnees$dtow,donnees$invertebrate,xlab='distance chalut�e',ylab='biomasse d invert�br�s') 
plot(donnees$depth, donnees$temperature,xlab='profondeur', ylab='temperature')
boxplot(donnees$invertebrate~donnees$strate)
boxplot(donnees$invertebrate~donnees$year) 
boxplot(donnees$invertebrate~donnees$sediment)
boxplot(donnees$invertebrate~donnees$starfish)
boxplot(donnees$invertebrate~donnees$urchin)

#ACP sur les variables quantitatives

donnees_ACP=donnees[,c(3,4,5,6,7,8,10)] #S�lection des variables quantitatives
ACP=PCA(donnees_ACP, scale.unit = TRUE, ncp = 7,)

#Regression lineaire multiple sur les variables quantitatives

reglm=lm(invertebrate~longitude+latitude+dtow+depth+temperature+totconsum,data=donnees_1) 
par(mfrow=c(2,2))
plot(reglm) #hypoth�ses du mod�le
summary(reglm) #Test du mod�le et estimateurs des variables
Cov=vcov(reglm) #Correlation entre les estimateurs de param�tre
Cor=cov2cor(Cov)
kable(Cor)
reg0=lm(invertebrate~1,data=donnees_1)
select.variables=step(reg0,scope=~latitude+longitude+dtow+depth+temperature+totconsum,direction ="both" )
summary(select.variables)

#Correction du mod�le:
regls=lm(invertebrate~dtow,data=donnees_1)#test de si la distance de chalutage est significative sur la quantit� p�ch�e
summary(regls)
par(mfrow=c(2,2))
plot(regls)
donnees=donnees[-c(326,327,410),] #retrait point aberrants
donnees_1=na.omit(donnees)
ratio=(donnees_1$invertebrate+1)/(donnees_1$dtow) #prise en compte de la distance de chalutage dans la quantit� d'invert�br�s captur�s

reglm_corr_log=lm(log(ratio)~longitude+latitude+depth+temperature+totconsum,data=donnees_1) #passage au logarithme du mod�le pour le corriger
par(mfrow=c(2,2))
plot(reglm_corr_log) #hypoth�ses du mod�le
summary(reglm_corr_log) #Test du mod�le et estimateurs des variables

Cov=vcov(reglm_corr_log) #Correlation entre les estimateurs de param�tre
Cor=cov2cor(Cov)
kable(Cor)
reg0=lm(invertebrate~1,data=donnees_1)
select.variables=step(reglm_corr_log )
summary(select.variables)

#Etude de l'effet du type de s�diments sur la quantit� d'invert�br�s:

reg_sediment=lm(log(ratio)~sediment,data=donnees_1) 
par(mfrow=c(2,2))
plot(reg_sediment)
summary(reg_sediment)

sediment.emmeans=emmeans(reg_sediment,pairwise ~ sediment,adjust="bonf")
sediment.emmeans

anova(reg_sediment)
Anova(reg_sediment)

#Etude de l'effet de la strate et de l'ann�e sur la quantit� d'invert�br�s

reg_strate_annee=lm(log(ratio)~strate*year,data=donnees_1)
par(mfrow=c(2,2))
plot(reg_strate_annee)
summary(reg_strate_annee)
anova(reg_strate_annee)
Anova(reg_strate_annee)

#Test d'un mod�le global avec Temp�rature, profondeur, strate et s�diments:

modele_global = lm(log(ratio)~depth+temperature+strate*sediment,data=donnees_1)
par(mfrow=c(2,2))
plot(modele_global)
summary(modele_global)
anova(modele_global)
Anova(modele_global)

#Application du mod�le � des esp�ces en particulier les oursins et les �toiles de mer:
test_oursin=glm(urchin~temperature+depth+totconsum+sediment,data=donnees_1, family = binomial(link = "logit"))
plot(test_oursin)
summary(test_oursin)
Anova(test_oursin)

test_etoile=glm(starfish~temperature+depth+totconsum+sediment,data=donnees_1, family = binomial(link = "logit"))
plot(test_etoile)
summary(test_etoile)
Anova(test_etoile)