## ----setup, include=FALSE-------------------------------------------------------------------------
fig_height <- 8
knitr::opts_chunk$set(echo = TRUE, fig.height = fig_height, 
                      eval = TRUE,
                      fig.width = 4/3 * fig_height,
                      comment = NA)


## ----librairires, message = FALSE-----------------------------------------------------------------
library(tidyverse)


## ----chargement_donnnees--------------------------------------------------------------------------
donnees <-read.csv(file="fraseriverdata.csv",
                   header = T, sep = ';', dec = ',')
knitr::kable(donnees)


## ----graphe_donnees-------------------------------------------------------------------------------
donnees %>% 
  select(year, spawners, total.run) %>% 
  rename(Spawners = spawners, Recrues = total.run) %>% 
  gather(-year, key = "Population", value = "Nombre") %>% 
  ggplot() + # Donnees représentées
  aes(x = year, y = Nombre, color = Population) +
  geom_line() +
  geom_point() +
  labs(x = "Année", y = "Nombre (milliers)", 
       title = "Evolution du nombre de saumons")


## ----graphe_stock_recrutement---------------------------------------------------------------------
nb_annees <- nrow(donnees) # Nombre d'années
stock_recrues <- tibble(Recrues = donnees$total.run[-1],
                        Stock = donnees$spawners[-nb_annees])
ggplot(stock_recrues) +
  aes(x = Stock, y = Recrues) +
  geom_point()


## ----fonctions_stock_recrutement------------------------------------------------------------------
ricker <- function(stock, alpha, beta){
  recrues <- alpha * stock * exp(-beta * stock)
  return(recrues)
}
beverton <- function(stock, alpha, beta){
  recrues <- alpha * stock / (beta + stock)
}


## ----parametres_exemple_stock_recrues, eval = T---------------------------------------------------
alpha_ricker <- 5; beta_ricker <- 0.15
alpha_beverton <- alpha_ricker; beta_beverton <- 5


## ----exemple_stock_recrutement--------------------------------------------------------------------

stock <- seq(0, 100, length.out = 1001)
tibble(Stock = stock) %>% 
  mutate(Ricker = ricker(Stock, alpha_ricker, beta_ricker),
         Beverton = beverton(Stock, alpha_beverton, beta_beverton)) %>% 
  gather(-Stock, key = "Modele", value = "Recrues") %>% 
  ggplot(aes(x = Stock, y = Recrues, col = Modele)) +
  geom_line()


## ----graphique_msy, warning = F-------------------------------------------------------------------
table_indicateurs <- tibble(Stock = c(1,
                                      0.5 * log(alpha_ricker) - 
                                        0.07 * log(alpha_ricker)^2) / beta_ricker) %>%
  mutate(Recrue = ricker(Stock, alpha_ricker, beta_ricker))
  
tibble(Stock = seq(0, 10, length.out = 101)) %>% 
  mutate(Recrue = ricker(Stock, alpha_ricker, beta_ricker)) %>% 
  ggplot(aes(x = Stock, y = Recrue)) +
  geom_line() +
  geom_abline(slope = 1, intercept = 0, linetype = 3) +
  geom_segment(data = table_indicateurs, 
               aes(xend = Stock), yend = 0, linetype = 2, col = "red") + 
  geom_segment(data = table_indicateurs, 
               aes(yend = Recrue), xend = 0, linetype = 2, col = "red") +
  annotate("text", 
           x = c(0, 0, table_indicateurs$Stock, 10), 
           y = c(table_indicateurs$Recrue, 0, 0, 10), 
           label = c(expression(R[max]), expression(R[MSY]),
                     expression(S[max]), expression(S[MSY]),
                     expression(R==S)), size = 8) 


## ----premiere_partie_solution, echo = T-----------------------------------------------------------
# On commence par bien importer les données
donnees_completes <-read.csv(file="fraseriverdata.csv",
                             header = T, sep = ';', dec = ',')
nb_annees <- nrow(donnees_completes) # Nombre d'années

# On créé un tableau stock recrutement supplémentaire
stock_recrue_vraies <- data.frame(Stock = donnees$spawners[-nb_annees], # On enlève la dernière année
                                   Recrue = donnees$total.run[-1], # On enlève la première année
                                  Annee = donnees$year[-1] # L'année est celle du recrutement. Le stock associé est celui de l'année n - 2
                        )


## ----premier_modele_jags--------------------------------------------------------------------------
script_modele <-"
model
{
	#sampling distribution (likelihood);
  for( i in 2 : N ) 					
  {
    Rm[i] <- log(alpha*Sobs[i-1]*exp(-beta*Sobs[i-1]));
    Robs[i ] ~ dlnorm(Rm[i],  isigmaR2);
    resids[i] <- Robs[i ] / exp(Rm[i]);
    lnresids[i] <- log(resids[i]);
  } # End loop
  
  #prior distribution of parameters;						
	
  alpha ~ dunif(0,30); # Loi non informative 
	beta ~ dunif(0,1); # Taux de décroissance
	isigmaR2 ~ dunif(0, 10); # Prior sur la précision
	sigmaR2 <- 1 / isigmaR2;

  #management quantities to be saved;				
	# Formules données dans le cours

  Smax <- 1 / beta;
	Rmax <- alpha / beta*exp(-1.0);
	SMSY <- log(alpha) / beta*(0.5 - 0.07 * log(alpha));
	uMSY <- 0.5 * log(alpha) - 0.07 * log(alpha) * log(alpha);
	MSY <- alpha * SMSY * exp(-beta * SMSY) - SMSY;
		
	} # End model
"


## ----burnnin_jags_modele1-------------------------------------------------------------------------
# Données
library(rjags) # Pour l'ajustement
donnees_jags <- list(N = nb_annees, # Nombre N
                     Robs = donnees$total.run, # Variable réponse
                     Sobs = donnees$spawners) # Variable explicative
n_chains <- 3 # Nombre de chaînes MCMC lancées
n_burn <- 5000 # Nombre d'itérations à oublier (on admettra qu'on a atteint
# la loi stationnaire ensuite)

# Points de départ des différentes chaînes
pt_depart1 <- list(alpha = 10, beta =0.5, isigmaR2 = 0.1)
pt_depart2 <- list(alpha = 15, beta =0.6, isigmaR2 = 0.4)
pt_depart3 <- list(alpha = 5, beta =0.4, isigmaR2 = 0.5)

# Initialisation du modèle
premier_modele <- jags.model(file=textConnection(script_modele), 
                     data = donnees_jags, 
                     inits = list(pt_depart1, pt_depart2, pt_depart3),
                     n.chains = n_chains, 
                     n.adapt = n_burn, # On lance n_burn iteration pour atteindre la loi stationnaire
                     quiet = FALSE) # Montre l'avancement


## ----echantillonnage_jags_modele1-----------------------------------------------------------------
n_iter <- 5000 # Nombre d'iterations effectuées
n_thin <- 10 # Pas entre deux itérations conservées
# Variables à conserver
variables_conservees <- c("alpha", "beta", "sigmaR2",
                         "MSY", "SMSY", "Rm", "lnresids")
echantillons_posterior_mcmc <- coda.samples(model = premier_modele, 
                                            variable.names = variables_conservees, 
                                            n.iter = n_iter,
                                            thin = n_thin)


## ----transformation_tibble------------------------------------------------------------------------
echantillons_posterior <- echantillons_posterior_mcmc %>% 
  map_dfr(# On applique aux 3 éléments de la liste la même function
    function(tab_mcmc){ # tab mcmc
      resultat <- tab_mcmc %>% # On prend le tableau mcmc
        as.data.frame() %>% # On transforme en data.frame
        mutate(num_echantillon = 1:nrow(tab_mcmc)) # On créé un indice de numéro
    },
    .id = "chaine") %>% # On créée une colonne chaine qui stocke le numero de chaine
  as_tibble()


## ----en_tete_tableau------------------------------------------------------------------------------
echantillons_posterior


## ----representation_chaine------------------------------------------------------------------------
echantillons_posterior %>% # Dans le tableau initial
  dplyr::select(num_echantillon, chaine, alpha, beta) %>%  # On sélectionne ces 4 colonnes
  gather(-num_echantillon, -chaine, key = "Parametre", # On transforme le tableau
         value = "valeur") %>% # en "tableau long" (regardez le résulat!)
  ggplot() + # On représente le résultat
  aes(x = num_echantillon, y = valeur, color = factor(chaine)) + # Les esthétiques
  geom_line() + # On trace 
  labs(x = "Numéro d'échantillon", y = "Valeur échantillonnée", # Habillage
       color = "Chaîne") +
  facet_wrap(~ Parametre, # Un graphe par paramètre
             scales = "free_y",
             labeller = label_parsed) # Pour que alpha soit transformé en lettre grecque


## ----gelman_plot----------------------------------------------------------------------------------
gelman.plot(echantillons_posterior_mcmc[, c("alpha","beta")])


## ----comparaison_prior_modele1--------------------------------------------------------------------
prior_alpha <- function(x){
  dunif(x, 0, 20) # Prior choisi dans le modèle
}
prior_beta <- function(x){
  dunif(x, 0, 1) # Prior choisi dans le modèle
}
graphe_alpha <- echantillons_posterior %>% # Dans le tableau initial
  dplyr::select(alpha) %>%  # On sélectionne ces 2 colonnes
  ggplot(aes(x = alpha, color = "Posterior")) +
  geom_density() + 
  labs(y = "Densité", x = "Valeur", title = expression(alpha),
       colour = "") +
  stat_function(mapping = aes(color = "Prior"),
                              fun = prior_alpha, xlim = c(0, 20))
graphe_beta <- echantillons_posterior %>% # Dans le tableau initial
  dplyr::select(beta) %>%  # On sélectionne ces 2 colonnes
  ggplot(aes(x = beta, color = "Posterior")) +
  geom_density()  +
  stat_function(fun = prior_alpha, xlim = c(0, 1), 
                mapping = aes(color = "Prior")) +
  labs(y = "Densité", x = "Valeur", title = expression(beta),
       colour = "")
gridExtra::grid.arrange(graphe_alpha, graphe_beta)


## ----visualisation2d------------------------------------------------------------------------------
p <- ggplot(echantillons_posterior, 
            aes(x = alpha, y = beta)) +
      geom_point(col = "red") + 
      geom_density_2d(col = "black") +
      theme(legend.position="none") +
  labs(x = expression(alpha), y = expression(beta))
ggExtra::ggMarginal(p, type = "histogram")


## ----echantillon_stock_recrue---------------------------------------------------------------------
echantillon_stock_recrue <- echantillons_posterior %>% 
  dplyr::select_at(vars("num_echantillon", "chaine",
                        starts_with("R"))) %>%  # Choix des variables de log residus
  gather(-num_echantillon, -chaine,
         key = "Annee", value = "Recrue", factor_key = T) %>% # On en fait un tableau long
  # (Regardez le résulat!)
  mutate(Annee = factor(Annee, labels = stock_recrue_vraies$Annee), # Transformation en années (facteur)
         Annee = as.numeric(as.character(Annee)), # Transformation en numérique
         Recrue = exp(Recrue)) %>% # Retour au monde réel
  left_join(y = stock_recrue_vraies[, c("Annee", "Stock")])


## ----graphique_predictions_recrues_stock----------------------------------------------------------
ggplot(echantillon_stock_recrue,
       aes(x = Stock, y = Recrue)) +
  geom_boxplot(aes(group = Stock,
                 color = "Echantillon posterior"), alpha = 0.1) +
  geom_line(data = stock_recrue_vraies, 
            mapping = aes(color = "Observations")) +
  labs(title = "Trajectoires du recrutement", y = "Recrues",
       color = "") +
  scale_color_manual(values = c("darkgreen", "red"))


## ----ic_90----------------------------------------------------------------------------------------
echantillon_stock_recrue %>% 
  group_by(Annee, Stock) %>% 
  summarise(Mediane = median(Recrue),
            Borne_sup = quantile(Recrue, prob = 0.95),
            Borne_inf = quantile(Recrue, prob = 0.05)) %>% 
  ggplot(aes(x = Stock)) +
  geom_ribbon(aes(ymin = Borne_inf, ymax = Borne_sup), 
              alpha=0.5, fill = "lightblue") +
  geom_line(aes(y = Mediane)) +
  geom_point(data = stock_recrue_vraies, # On représente les vraies données 
             mapping =  aes(y = Recrue)) +
  labs(y = "Recrue", title = "Intervalle de confiance à 90%")


## ----graphique_predictions_recrues_temps----------------------------------------------------------
echantillons_posterior %>% 
  dplyr::select_at(vars("num_echantillon", "chaine",
                        starts_with("R"))) %>%  # Choix des variables de log residus
  gather(-num_echantillon, -chaine,
         key = "Annee", value = "Valeur", factor_key = T) %>% # On en fait un tableau long
  # (Regardez le résulat!)
  mutate(Annee = factor(Annee, labels = donnees$year[-1]),
         Annee = as.numeric(as.character(Annee)),
         Valeur = exp(Valeur)) %>% 
  ggplot(aes(x = Annee, y = Valeur)) + 
  geom_line(aes(group = interaction(num_echantillon, chaine),
                color = "Echantillon posterior"), alpha = 0.1) +
  geom_line(data = donnees, mapping = aes(x = year, y = total.run, color = "Observations")) +
  labs(title = "Trajectoires du recrutement", y = "Recrues",
       color = "") +
  scale_color_manual(values = c("darkgreen", "red"))


## ----graphique_log_residus------------------------------------------------------------------------
echantillons_posterior %>% 
  dplyr::select_at(vars(starts_with("lnresids"))) %>%  # Choix des variables de log residus
  gather(key = "Annee", value = "Valeur", factor_key = T) %>% # On en fait un tableau long
  # (Regardez le résulat!)
  mutate(Annee = factor(Annee, labels = donnees$year[-1])) %>% 
  ggplot(aes(x = Annee, y = Valeur)) + 
  geom_boxplot() +
  labs(title = "Distribution du log des résidus")


## ----script_modele_hmm, echo = T------------------------------------------------------------------
script_second_modele <- "
model{
  #### INITIALISATION
  
  ## Modele de DYNAMIQUE 
  # Premieres valeurs de mortalite par peche et recrutement
	FO	 ~ dlnorm(1.3, 10); # centree sur 4, de variance 1.25
	RO	 ~ dlnorm(2, 10); # Centree sur 8.1, de variance 3
  Rm[1] <- log(alpha*exp(-FO)*RO*exp(-beta*exp(-FO)*RO));
  Fm[1] <- log(FO);
  F[1] ~ dlnorm(Fm[1], isigmaF2);
  R[1] ~ dlnorm(Rm[1], isigmaR2);
  C[1] <- log(R[1] * (1-exp(-F[1])));
	S[1] <- log(R[1] * exp(-F[1]));
  
  ## Modele d'OBSERVATION
  Cobs[1] ~ dlnorm(C[1], itauC2);
	Sobs[1] ~ dlnorm(S[1], itauS2);
  
  ## Résidus
  residsS[1] <- Sobs[1] - exp(S[1]);
  residsC[1] <- Cobs[1] - exp(C[1]);
  
  #### PROPAGATION
  for( i in 2 : N )
  {
  ## Modele de DYNAMIQUE 
    Fm[i] <- log(F[i-1]); # Moyenne de la nouvelle mortalite par peche
    F[i ] ~ dlnorm(Fm[i],isigmaF2); # Nouvelle mortalité
    # Nouvelle moyenne de recrutement
    Rm[i] <- log(alpha*exp(-F[i-1])*R[i-1]*exp(-beta*exp(-F[i-1])*R[i-1]));
    R[i ] ~ dlnorm(Rm[i],isigmaR2); # Nouveau recrutement
    C[i] <- log(R[i]*(1-exp(-F[i]))); # On actualise les captures
    S[i]<- log(R[i]*exp(-F[i])); # On actualise le stock

  ## Modele d'OBSERVATION
    Cobs[i ] ~ dlnorm(C[i],itauC2); # Captures observees
    Sobs[i ] ~ dlnorm(S[i],itauS2); # Reproducteurs observes

  ## Résidus
    residsS[i] <- Sobs[i ]-exp(S[i]); # Pour S
    residsC[i] <- Cobs[i ]-exp(C[i]); # Pour C
    resids[i] <- (R[i])/exp(Rm[i]); # Pour R
    lnresids[i] <- log(resids[i]); # Pour log R
	} # Fin de la boucle sur le temps

  ## Paramètres des priors
	alpha ~ dunif(0, 20); # On informe un peu plus que précédemment
	beta ~ dunif(0,1);
		
	sigmaR2 ~ dunif(0, 1);
	isigmaR2 <- 1 / sigmaR2;
	sigmaF2 ~ dunif(0, 0.1); # Pas trop variable
	isigmaF2 <- 1 / sigmaF2;
  tauC2 ~ dunif(0, 1);
	itauC2 <- 1 / tauC2;
  tauS2 ~ dunif(0, 1);
	itauS2 <- 1 / tauS2;
		
  #management quantities;					
	Smax <- 1 / beta;
	Rmax <- alpha / beta*exp(-1.0);
	SMSY <- log(alpha)/ beta*(0.5-0.07*log(alpha));
	uMSY <- 0.5*log(alpha)-0.07*log(alpha)*log(alpha);
	MSY <- alpha*SMSY*exp(-beta*SMSY)-SMSY;
		
}
"


## ----burnnin_jags_modele2, echo = T---------------------------------------------------------------
# Données
donnees_modele2 <- list(N = nb_annees, # Nombre N
                     Cobs = donnees$catch, # Variable réponse
                     Sobs = donnees$spawners) # Variable explicative
n_chains <- 3 # Nombre de chaînes MCMC lancées
n_burn <- 10000 # Nombre d'itérations à oublier (on admettra qu'on a atteint
# la loi stationnaire ensuite)

# Initialisation du modèle
second_modele <- jags.model(file=textConnection(script_second_modele), 
                     data = donnees_modele2, 
                     n.chains = n_chains, 
                     n.adapt = n_burn, # On lance n_burn iteration pour atteindre la loi stationnaire
                     quiet = FALSE) # Montre l'avancement


## ----echantillonnage_jags_modele2, echo = T-------------------------------------------------------
n_iter <- 10000 # Nombre d'iterations effectuées
n_thin <- 10 # Pas entre deux itérations conservées

# Variables à conserver

variables_conservees <- c("alpha","beta", # Parametres de dynamique
                          "sigmaR2", "sigmaF2", "tauC2", "tauS2", # Parametres de variance
                          "R","F", # Dynamiques cachees
                          "MSY","SMSY", # Paramètres de gestion
                          "lnresids") # Residus
echantillons_posterior_mcmc_modele2 <- coda.samples(model = second_modele, 
                                            variable.names = variables_conservees, 
                                            n.iter = n_iter,
                                            thin = n_thin)


## ----echantillons_posterior_m2--------------------------------------------------------------------
echantillons_posterior_modele2 <- echantillons_posterior_mcmc_modele2 %>% 
  map_dfr(# On applique aux 3 éléments de la liste la même function
    function(tab_mcmc){ # tab mcmc
      resultat <- tab_mcmc %>% # On prend le tableau mcmc
        as.data.frame() %>% # On transforme en data.frame
        mutate(num_echantillon = 1:nrow(tab_mcmc)) # On créé un indice de numéro
    },
    .id = "chaine") %>% # On créée une colonne chaine qui stocke le numero de chaine
  as_tibble()

