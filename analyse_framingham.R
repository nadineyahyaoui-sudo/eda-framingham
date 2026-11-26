# =============================================================================
# Projet : Analyse Exploratoire de Données & Data Storytelling
# Données : Framingham Heart Study (maladie coronarienne à 10 ans)
# Problématique : Quels facteurs (âge, tension, cholestérol, tabac, glycémie)
#                 sont associés au risque de maladie coronarienne à 10 ans ?
# Utilisation : ouvrir ce fichier dans RStudio puis cliquer sur "Source"
#               (les figures sont enregistrées dans le dossier "figures/")
# =============================================================================


# 1. PACKAGES ----------------------------------------------------------------
# On installe automatiquement les packages manquants (script auto-porteur)
paquets <- c("tidyverse", "scales")
manquants <- paquets[!paquets %in% rownames(installed.packages())]
if (length(manquants) > 0) install.packages(manquants)

library(tidyverse)   # dplyr, tidyr, ggplot2, readr...
library(scales)      # affichage des pourcentages

dir.create("figures", showWarnings = FALSE)   # dossier pour les graphiques


# 2. IMPORT DES DONNÉES ------------------------------------------------------
url <- paste0(
  "https://raw.githubusercontent.com/RiahiVahid/",
  "Framingham-Heart-study-dataset/master/framingham_heart_disease.csv"
)
df_brut <- read_csv(url, show_col_types = FALSE)

cat("Dimensions :", nrow(df_brut), "lignes x", ncol(df_brut), "colonnes\n")
glimpse(df_brut)


# 3. DONNÉES MANQUANTES ------------------------------------------------------
# Nombre et pourcentage de NA par variable
na_resume <- df_brut %>%
  summarise(across(everything(), ~ sum(is.na(.x)))) %>%
  pivot_longer(everything(), names_to = "variable", values_to = "nb_na") %>%
  mutate(pct_na = round(100 * nb_na / nrow(df_brut), 1)) %>%
  filter(nb_na > 0) %>%
  arrange(desc(nb_na))

print(na_resume)

# Les NA représentent peu de lignes : on les retire seulement dans les
# analyses qui en ont besoin (drop_na ciblé), on ne supprime rien globalement.


# 4. PRÉPARATION DES VARIABLES -----------------------------------------------
# Variables catégorielles avec des libellés lisibles, et variables numériques
# rééchelonnées pour interpréter les odds ratios "par 10 unités".
df <- df_brut %>%
  mutate(
    chd    = factor(TenYearCHD, levels = c(0, 1),
                    labels = c("Sans maladie", "Maladie")),
    fumeur = factor(currentSmoker, levels = c(0, 1),
                    labels = c("Non-fumeur", "Fumeur")),
    sexe   = factor(male, levels = c(0, 1), labels = c("Femme", "Homme")),
    age10     = age / 10,       # effet par 10 ans
    sysBP10   = sysBP / 10,     # effet par 10 mmHg
    totChol10 = totChol / 10,   # effet par 10 mg/dL
    glucose10 = glucose / 10    # effet par 10 mg/dL
  )


# 5. VALEURS ABERRANTES (règle de l'IQR) -------------------------------------
# Une valeur est signalée si elle sort de [Q1 - 1.5*IQR ; Q3 + 1.5*IQR].
# On les SIGNALE mais on les garde : ce sont des valeurs médicalement possibles
# (ex. glycémie 386 chez un diabétique). On utilise donc des statistiques
# robustes (médiane, test de Wilcoxon).
signale_iqr <- function(x) {
  q   <- quantile(x, c(0.25, 0.75), na.rm = TRUE)
  iqr <- q[2] - q[1]
  x < q[1] - 1.5 * iqr | x > q[2] + 1.5 * iqr
}

aberrantes <- df %>%
  summarise(across(c(age, sysBP, diaBP, totChol, BMI, glucose),
                   ~ sum(signale_iqr(.x), na.rm = TRUE))) %>%
  pivot_longer(everything(), names_to = "variable", values_to = "nb_aberrantes")

print(aberrantes)


# 6. STATISTIQUES DESCRIPTIVES PAR GROUPE ------------------------------------
tableau_desc <- df %>%
  group_by(chd) %>%
  summarise(across(
    c(age, sysBP, totChol, BMI, glucose),
    list(
      moy = ~ mean(.x, na.rm = TRUE),
      med = ~ median(.x, na.rm = TRUE),
      sd  = ~ sd(.x, na.rm = TRUE),
      q1  = ~ quantile(.x, 0.25, na.rm = TRUE),
      q3  = ~ quantile(.x, 0.75, na.rm = TRUE)
    ),
    .names = "{.col}_{.fn}"
  )) %>%
  pivot_longer(-chd, names_to = c("variable", "stat"), names_sep = "_") %>%
  pivot_wider(names_from = stat, values_from = value) %>%
  mutate(across(where(is.numeric), ~ round(.x, 1)))

print(tableau_desc, n = Inf)
write_csv(tableau_desc, "figures/tableau_descriptif.csv")

# Prévalence globale de la maladie
cat("\nTaux de maladie à 10 ans :",
    percent(mean(df$TenYearCHD), accuracy = 0.1), "\n")


# 7. GRAPHIQUES --------------------------------------------------------------
# Thème commun, grande police pour la lisibilité sur poster A0
theme_poster <- theme_minimal(base_size = 20) +
  theme(plot.title = element_text(face = "bold"),
        legend.position = "bottom")

couleurs <- c("Sans maladie" = "#2A9D8F", "Maladie" = "#E63946")

# Graphique 1 : distribution de l'âge selon la maladie (histogramme)
g1 <- ggplot(df, aes(x = age, fill = chd)) +
  geom_histogram(binwidth = 2, position = "identity", alpha = 0.6,
                 colour = "white") +
  scale_fill_manual(values = couleurs) +
  labs(title = "Les malades sont plus âgés",
       x = "Âge (années)", y = "Nombre de patients", fill = NULL) +
  theme_poster

# Graphique 2 : tension systolique selon la maladie (boxplot)
g2 <- ggplot(df, aes(x = chd, y = sysBP, fill = chd)) +
  geom_boxplot(alpha = 0.7, outlier.alpha = 0.4) +
  scale_fill_manual(values = couleurs) +
  labs(title = "Tension systolique plus élevée chez les malades",
       x = NULL, y = "Tension systolique (mmHg)") +
  theme_poster + theme(legend.position = "none")

# Graphique 3 : taux de maladie selon le tabagisme (barplot)
taux_tabac <- df %>%
  group_by(fumeur) %>%
  summarise(taux = mean(TenYearCHD), n = n())

g3 <- ggplot(taux_tabac, aes(x = fumeur, y = taux, fill = fumeur)) +
  geom_col(width = 0.6) +
  geom_text(aes(label = percent(taux, accuracy = 0.1)),
            vjust = -0.5, size = 8) +
  scale_y_continuous(labels = percent, limits = c(0, 0.25)) +
  scale_fill_manual(values = c("Non-fumeur" = "#457B9D",
                               "Fumeur" = "#E76F51")) +
  labs(title = "Taux de maladie selon le tabagisme",
       x = NULL, y = "Taux de maladie à 10 ans") +
  theme_poster + theme(legend.position = "none")

# Graphique 4 : âge vs tension, avec tendance (nuage de points)
g4 <- ggplot(df, aes(x = age, y = sysBP, colour = chd)) +
  geom_point(alpha = 0.35, size = 1.5) +
  geom_smooth(method = "lm", se = TRUE, linewidth = 1.2) +
  scale_colour_manual(values = couleurs) +
  labs(title = "Tension et âge",
       x = "Âge (années)", y = "Tension systolique (mmHg)", colour = NULL) +
  theme_poster

# Sauvegarde en haute résolution (utilisables dans le poster LaTeX)
ggsave("figures/g1_age_histogramme.png",  g1, width = 10, height = 7, dpi = 300)
ggsave("figures/g2_tension_boxplot.png",  g2, width = 10, height = 7, dpi = 300)
ggsave("figures/g3_tabac_barplot.png",    g3, width = 10, height = 7, dpi = 300)
ggsave("figures/g4_age_tension_nuage.png", g4, width = 10, height = 7, dpi = 300)


# 8. TESTS STATISTIQUES ------------------------------------------------------
# Test 1 : la tension systolique diffère-t-elle entre malades et non-malades ?
# Wilcoxon (non paramétrique) car la distribution est asymétrique et il y a
# des valeurs extrêmes.
# H0 : même distribution de sysBP dans les deux groupes.
test_wilcoxon <- wilcox.test(sysBP ~ chd, data = df)
print(test_wilcoxon)

# Test 2 : le tabagisme est-il lié à la maladie ? (chi-deux d'indépendance)
# H0 : tabagisme et maladie sont indépendants.
table_tabac <- table(df$fumeur, df$chd)
print(table_tabac)
test_chi2 <- chisq.test(table_tabac)
print(test_chi2)


# 9. RÉGRESSION LOGISTIQUE ---------------------------------------------------
# Modèle : P(maladie) en fonction de l'âge, la tension, le cholestérol, le
# tabac et la glycémie.
#   logit(p) = b0 + b1*age10 + b2*sysBP10 + b3*totChol10 + b4*fumeur + b5*glucose10
# On ne garde que les lignes complètes pour ces variables.
df_modele <- df %>%
  select(TenYearCHD, age10, sysBP10, totChol10, currentSmoker, glucose10) %>%
  drop_na()

cat("\nObservations utilisées dans le modèle :", nrow(df_modele),
    "sur", nrow(df), "\n")

modele <- glm(TenYearCHD ~ age10 + sysBP10 + totChol10 +
                currentSmoker + glucose10,
              data = df_modele, family = binomial)

print(summary(modele))

# Odds ratios et intervalles de confiance à 95 %
# OR > 1 : le facteur augmente la probabilité de maladie (toutes choses égales)
resultats_or <- tibble(
  variable = names(coef(modele)),
  OR       = exp(coef(modele)),
  IC_bas   = exp(coef(modele) - 1.96 * sqrt(diag(vcov(modele)))),
  IC_haut  = exp(coef(modele) + 1.96 * sqrt(diag(vcov(modele)))),
  p_value  = summary(modele)$coefficients[, 4]
) %>%
  filter(variable != "(Intercept)") %>%
  mutate(across(where(is.numeric), ~ round(.x, 3)))

print(resultats_or)
write_csv(resultats_or, "figures/odds_ratios.csv")

# Qualité du modèle : pseudo-R² de McFadden
modele_nul <- glm(TenYearCHD ~ 1, data = df_modele, family = binomial)
r2_mcfadden <- 1 - as.numeric(logLik(modele) / logLik(modele_nul))
cat("Pseudo-R2 de McFadden :", round(r2_mcfadden, 3), "\n")

# Graphique 5 : forest plot des odds ratios
noms_clairs <- c(age10 = "Âge (+10 ans)",
                 sysBP10 = "Tension (+10 mmHg)",
                 totChol10 = "Cholestérol (+10 mg/dL)",
                 currentSmoker = "Fumeur (vs non)",
                 glucose10 = "Glycémie (+10 mg/dL)")

g5 <- resultats_or %>%
  mutate(variable = recode(variable, !!!noms_clairs)) %>%
  ggplot(aes(x = OR, y = reorder(variable, OR))) +
  geom_vline(xintercept = 1, linetype = "dashed") +
  geom_pointrange(aes(xmin = IC_bas, xmax = IC_haut),
                  linewidth = 1.2, size = 1, colour = "#E63946") +
  labs(title = "Odds ratios (IC à 95 %)",
       x = "Odds ratio", y = NULL) +
  theme_poster

ggsave("figures/g5_odds_ratios.png", g5, width = 10, height = 7, dpi = 300)

# Rappel d'interprétation : étude d'observation, donc on parle d'ASSOCIATION,
# jamais de causalité.
cat("\nAnalyse terminée. Figures et tableaux dans le dossier 'figures/'.\n")