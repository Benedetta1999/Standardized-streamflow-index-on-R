pacman::p_load(readr,SPEI,openxlsx,tidyverse,fitdistrplus,actuar,evd,dyplr,lubridate,gridExtra)

setwd("C:\\Users\\Benedetta\\OneDrive - Politecnico di Torino\\Documenti\\SSI\\Andonno Gesso")

my_df <- read_csv("ANDONNO GESSO.csv")

my_df$DATA <- as.Date(my_df$DATA, format="%Y-%m-%d")

my_df$Anno <- format(my_df$DATA, "%Y")
my_df$Mese <- format(my_df$DATA, "%m")

dati_portate_medie <- data.frame(Year = integer(), Month = integer(), prcp = numeric(), stringsAsFactors = FALSE)

####ELIMINANDO I MESI CON MENO DI 4 VALORI
# # Ciclo per calcolare le precipitazioni cumulative per ogni mese e anno
# for (anno in unique(my_df$Anno)) {
#   for (mese in 1:12) {
#     # Seleziona i dati per il mese e l'anno correnti
#     dati_mese_anno <- my_df[my_df$Mese == sprintf("%02d", mese) & my_df$Anno == anno, ]
#     
#     # Controlla se ci sono 4 o più valori mancanti nel mese corrente
#     if (sum(is.na(dati_mese_anno$"Portata fiume (m³/s)")) >= 4) {
#       media_mese <- NA
#     } else {
#       # Sostituisci i valori mancanti con la media del mese
#       media_mese <- mean(dati_mese_anno$"Portata fiume (m³/s)", na.rm = TRUE)
#       dati_mese_anno$"Portata fiume (m³/s)"[is.na(dati_mese_anno$"Portata fiume (m³/s)")] <- media_mese
#       
#       # Calcola la media della portata fiume per il mese corrente
#       media_mese <- mean(dati_mese_anno$"Portata fiume (m³/s)", na.rm = TRUE)
#     }
#     
#     # Aggiungi i dati al dataframe delle medie
#     dati_portate_medie <- rbind(dati_portate_medie, data.frame(Year = anno, Month = mese, portata_media = media_mese))
#   }
# }
# # Carica il pacchetto necessario per la manipolazione dei dati

####SENZA ELIMINARE I MESI CON MENO DI 4 VALORI

for (anno in unique(my_df$Anno)) {
  for (mese in 1:12) {
    
    dati_mese_anno <- my_df[my_df$Mese == sprintf("%02d", mese) & my_df$Anno == anno, ]
    
    
    media_mese <- mean(dati_mese_anno$"Portata fiume (m³/s)", na.rm = TRUE)
    dati_mese_anno$"Portata fiume (m³/s)"[is.na(dati_mese_anno$"Portata fiume (m³/s)")] <- media_mese
    
    
    media_mese <- mean(dati_mese_anno$"Portata fiume (m³/s)", na.rm = TRUE)
    
    
    dati_portate_medie <- rbind(dati_portate_medie, data.frame(Year = anno, Month = mese, portata_media = media_mese))
  }
}

dati_portate_medie$volume_medio <- dati_portate_medie$portata_media * 2592000

media_volume_medio <- mean(dati_portate_medie$volume_medio, na.rm = TRUE)


dati_portate_medie$media_volume_medio <- media_volume_medio

dati_portate_medie$Year <- as.numeric(dati_portate_medie$Year)

dati_portate_medie$Year_Month <- ymd(paste(dati_portate_medie$Year, dati_portate_medie$Month, "01", sep = "-"))

# Crea il grafico delle portate medie mensili
ggplot(dati_portate_medie, aes(x = Year_Month, y = portata_media)) +
  geom_line(group = 1, color = "black") +
  labs(title = "Portate Medie Mensili a GRAMO",
       x = "",
       y = "Portata Media (m³/s)") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  scale_x_date(date_breaks = "1 years", date_labels = "%Y")


# Creiamo un dataframe per le portate mensili
portate_mensili <- data.frame(Year = unique(dati_portate_medie$Year))  # Creiamo una colonna per gli anni

# Estraiamo i dati delle portate medie per ogni mese
for (mese in 1:12) {
  col_name <- as.character(mese)  # Nome della colonna per il mese
  portate_mensili[[col_name]] <- rep(NA, nrow(portate_mensili))  # Inizializziamo la colonna con NA
  
  for (anno in unique(dati_portate_medie$Year)) {
    # Troviamo il valore della portata media per il mese e l'anno corrente
    valore_portata <- dati_portate_medie$portata_media[dati_portate_medie$Year == anno & dati_portate_medie$Month == mese]
    
    # Verifica se il vettore valore_portata è vuoto
    if (length(valore_portata) > 0) {
      # Assegniamo il valore alla colonna corrispondente solo se il vettore non è vuoto
      portate_mensili[[col_name]][portate_mensili$Year == anno] <- valore_portata
    }
  }
}

# Imposta gli anni come indice del dataframe portate_mensili
rownames(portate_mensili) <- portate_mensili$Year

# Inizializzazione dei fit per le distribuzioni
fit_ln <- list()
fit_weibull <- list()
fit_ll <- list()
fit_gev<-list()
fit_pp<-list()

# Loop attraverso ogni colonna di portate_mensili (eccetto la colonna dell'anno)
for (col in names(portate_mensili)[-1]) {
  # Filtra i valori non NA nella colonna corrente
  data <- portate_mensili[[col]][complete.cases(portate_mensili[[col]])]
  
  # Fit delle distribuzioni ai dati
  fit_ln[[col]] <- suppressWarnings(fitdist(data, "lnorm"))
  fit_weibull[[col]] <- suppressWarnings(fitdist(data, "weibull"))
  fit_ll[[col]] <- suppressWarnings(fitdist(data, "llogis"))
}

# Creazione della lista per le distribuzioni migliori per ogni mese
migliori_distribuzioni <- list()

# Loop attraverso ogni colonna di portate_mensili
for (col in names(portate_mensili)[-1]) {
  # Inizializzazione della lista per questo mese
  migliori_distribuzioni[[col]] <- list()
  
  # Estrai i valori AIC per ogni distribuzione
  aic_ln <- fit_ln[[col]]$aic
  aic_weibull <- fit_weibull[[col]]$aic
  aic_ll <- fit_ll[[col]]$aic
  
  # Trova la distribuzione con il valore AIC più basso
  distribuzione_migliore <- which.min(c(aic_ln, aic_weibull, aic_ll))
  
  # Assegna la distribuzione migliore a questo mese
  migliori_distribuzioni[[col]]$nome <- switch(distribuzione_migliore,
                                               "1" = "lnorm",
                                               "2" = "weibull",
                                               "3" = "llogis")
  
  # Memorizza i parametri della distribuzione migliore
  parametri_migliori <- switch(distribuzione_migliore,
                               "1" = fit_ln[[col]]$estimate,
                               "2" = fit_weibull[[col]]$estimate,
                               "3" = fit_ll[[col]]$estimate)
  
  # Aggiungi i parametri alla lista
  migliori_distribuzioni[[col]]$parametri <- parametri_migliori
}



# Inizializzazione della lista per le probabilità cumulative
prob_cumulative <- list()

# Loop attraverso ogni colonna di portate_mensili
for (col in names(portate_mensili)[-1]) {
  # Estrai la distribuzione migliore per questo mese
  distribuzione_migliore <- migliori_distribuzioni[[col]]$nome
  parametri <- migliori_distribuzioni[[col]]$parametri
  
  # Calcola la probabilità cumulata con la distribuzione migliore
  if (distribuzione_migliore == "lnorm") {
    prob_cumulative[[col]] <- plnorm(portate_mensili[[col]], meanlog = parametri[1], sdlog = parametri[2])
  } else if (distribuzione_migliore == "weibull") {
    prob_cumulative[[col]] <- pweibull(portate_mensili[[col]], shape = parametri[1], scale = parametri[2])
  } else if (distribuzione_migliore == "llogis") {
    prob_cumulative[[col]] <- pllogis(portate_mensili[[col]], shape = parametri[1], scale = parametri[2])
  }
}

# Visualizza la lista delle probabilità cumulative
print(prob_cumulative)


#######PLOT CHECK EMPIRICAL - THEORETICAL DISTRIBUTION FUNCTION

plot_cdf <- function(data, distribuzione, parametri, titolo) {
  empirical_cdf <- ecdf(data)
  x_values <- seq(min(data, na.rm = TRUE), max(data, na.rm = TRUE), length.out = 100)
  
  theoretical_cdf <- switch(distribuzione,
                            "lnorm" = plnorm(x_values, meanlog = parametri[1], sdlog = parametri[2]),
                            "weibull" = pweibull(x_values, shape = parametri[1], scale = parametri[2]),
                            "llogis" = pllogis(x_values, shape = parametri[1], scale = parametri[2]))
  
  ggplot() +
    stat_ecdf(data = data.frame(x = data), aes(x = x), geom = "point", color = "blue", linetype = "dashed", size = 1) +
    geom_line(data = data.frame(x = x_values, y = theoretical_cdf), aes(x = x, y = y), color = "red", size = 1) +
    labs(title = titolo, x = "Valori", y = "Probabilità cumulativa") +
    theme_minimal() +
    scale_color_manual(name = "Legenda", values = c("Densità empirica" = "blue", "Densità teorica" = "red")) +
    guides(color = guide_legend(override.aes = list(linetype = c("dashed", "solid"))))
}

# Crea una lista di tutti i plot
plots <- list()
for (col in names(portate_mensili)[-1]) {
  data <- portate_mensili[[col]][complete.cases(portate_mensili[[col]])]
  distribuzione <- migliori_distribuzioni[[col]]$nome
  parametri <- migliori_distribuzioni[[col]]$parametri
  titolo <- paste("CDF per", col, "con distribuzione", distribuzione)
  plots[[col]] <- plot_cdf(data, distribuzione, parametri, titolo)
}

# Visualizza tutti i plot in una griglia di 3 colonne
do.call(grid.arrange, c(plots, ncol = 3))



### CALCOLO PROB CUM NORM

# Inizializzazione della lista per le probabilità cumulative normalizzate
prob_cumulative_normalizzate <- list()

# Loop attraverso ogni colonna di portate_mensili
for (col in names(portate_mensili)[-1]) {
  # Estrai la probabilità cumulata per questo mese
  prob_cumulative_mese <- prob_cumulative[[col]]
  
  # Normalizza la probabilità cumulata con qnorm
  prob_cumulative_normalizzata <- qnorm(prob_cumulative_mese)
  
  # Assegna la probabilità cumulata normalizzata alla lista
  prob_cumulative_normalizzate[[col]] <- prob_cumulative_normalizzata
}

# Visualizza la lista delle probabilità cumulative normalizzate
print(prob_cumulative_normalizzate)



# Inizializzazione del dataframe vuoto con gli anni come indici delle righe
df_prob_cumulative <- data.frame(Year = as.integer(rownames(portate_mensili)), row.names = NULL, stringsAsFactors = FALSE)

# Unione dei dataframe della lista prob_cumulative_normalizzate
df_prob_cumulative <- cbind(df_prob_cumulative, do.call(cbind, prob_cumulative_normalizzate))

# Rinomina le colonne con i numeri dei mesi
colnames(df_prob_cumulative)[-1] <- 1:12

# Visualizza il dataframe risultante
print(df_prob_cumulative)

# Creazione del dataframe distribuzione_ogni_mese
distribuzione_ogni_mese <- data.frame(Year = integer(),
                                      Month = integer(),
                                      Portata_Media = numeric(),
                                      SSI = numeric(),
                                      stringsAsFactors = FALSE)

# Loop attraverso ogni riga del dataframe portate_mensili (ogni anno)
for (year in rownames(portate_mensili)) {
  # Loop attraverso ogni colonna del dataframe portate_mensili (ogni mese)
  for (month in names(portate_mensili)[-1]) {
    # Estrai il valore della portata media per questo mese e anno
    portata_media <- portate_mensili[year, month]
    
    # Estrai l'indice SSI corrispondente all'anno e mese corrente
    SSI_index <- as.integer(year) - min(df_prob_cumulative$Year) + 1
    
    # Estrai il valore di SSI corrispondente all'indice calcolato e al mese corrente
    SSI <- df_prob_cumulative[SSI_index, as.character(month)]
    
    # Aggiungi i dati al dataframe distribuzione_ogni_mese
    distribuzione_ogni_mese <- rbind(distribuzione_ogni_mese, data.frame(
      Year = as.integer(year),
      Month = as.integer(month),
      Portata_Media = portata_media,
      SSI = SSI
    ))
  }
}

# Visualizza il dataframe distribuzione_ogni_mese
print(distribuzione_ogni_mese)

write.xlsx(distribuzione_ogni_mese, "SSI_VinadioDemonte.xlsx")

###########EVENTI CON ONSET


# Inizializzazione della colonna degli eventi di siccità nel dataframe distribuzione_ogni_mese
distribuzione_ogni_mese$Evento_siccita_ssi <- ""

# Trova gli eventi di siccità SSI
for (i in 1:nrow(distribuzione_ogni_mese)) {
  if (!is.na(distribuzione_ogni_mese$SSI[i]) && distribuzione_ogni_mese$SSI[i] < -1) {
    # Controllo dei mesi precedenti
    for (j in i:1) {
      if (is.na(distribuzione_ogni_mese$SSI[j])) {
        break
      } else if (distribuzione_ogni_mese$SSI[j] < 0) {
        distribuzione_ogni_mese$Evento_siccita_ssi[j:i] <- "Siccita in corso SSI"
      } else {
        break
      }
    }
    # Controllo dei mesi successivi
    for (k in i:nrow(distribuzione_ogni_mese)) {
      if (is.na(distribuzione_ogni_mese$SSI[k])) {
        break
      } else if (distribuzione_ogni_mese$SSI[k] < 0) {
        distribuzione_ogni_mese$Evento_siccita_ssi[i:k] <- "Siccita in corso SSI"
      } else {
        break
      }
    }
  }
}

distribuzione_ogni_mese$Colore_ssi <- ifelse(distribuzione_ogni_mese$Evento_siccita_ssi == "Siccita in corso SSI", "Rosso", "Nero")

distribuzione_ogni_mese$Date <- as.Date(with(distribuzione_ogni_mese, paste(Year, Month, "01", sep = "-")), "%Y-%m-%d")

min_date <- min(distribuzione_ogni_mese$Date, na.rm = TRUE)
max_date <- max(distribuzione_ogni_mese$Date, na.rm = TRUE)

# Grafico dell'indice SSI con colori per gli eventi di siccità
grafico_ssi <- ggplot(distribuzione_ogni_mese, aes(x = Date, y = SSI, fill = Colore_ssi)) +
  geom_bar(stat = "identity", width = 25) + # Imposta una larghezza adeguata per le barre
  scale_fill_manual(values = c("Rosso" = "red", "Nero" = "black"), labels = c("Rosso" = "Drought Runs", "Nero" = "SSI")) +
  geom_hline(yintercept = c(-1), linetype = "dashed", color = "red", size = 0.5) +
  geom_hline(yintercept = c(1), linetype = "dashed", color = "blue", size = 0.5) +
  labs(x = "", y = "", title = "SSI indice a GRAMO") +
  theme_minimal() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    legend.position = c(0.93, 0.95),
    legend.title = element_blank(),  # Rimuove il titolo della legenda
    panel.grid.major = element_blank(),  # Rimuove la griglia principale
    panel.grid.minor = element_blank(),  # Rimuove la griglia secondaria
    axis.line = element_line(color = "black"),
    legend.box.background = element_rect(color = "black", size = 0.5),  # Aggiunge il riquadro attorno alla legenda
    panel.background = element_rect(fill = "white", color = "black"),  # Sfondo bianco per il pannello del grafico
    plot.background = element_rect(fill = "white", color = "black")    # Sfondo bianco per l'intera area del grafico
  ) +
  scale_x_date(
    expand = c(0, 0),
    limits = c(min_date, max_date),
    date_breaks = "1 years",
    date_labels = "%Y"
  )


print(grafico_ssi)

file_path <- "C:/Users/Benedetta/OneDrive - Politecnico di Torino/Documenti/SSI/Grafici SSI/SDEVI.png"

# Salva il grafico nel percorso specificato
ggsave(file_path, plot = grafico_ssi, width = 10, height = 6, dpi = 300)

# Crea una nuova dataframe per gli eventi siccitosi
eventi_siccitosi <- data.frame(
  Data_di_Inizio = character(),
  Data_di_Fine = character(),
  Durata = numeric(),
  stringsAsFactors = FALSE
)

# Trova gli eventi di siccità SSI
evento_in_corso <- FALSE
inizio_evento <- ""
fine_evento <- ""
durata_evento <- 0

for (i in 1:nrow(distribuzione_ogni_mese)) {
  if (distribuzione_ogni_mese$Colore_ssi[i] == "Rosso" && !evento_in_corso) {
    evento_in_corso <- TRUE
    inizio_evento <- paste(distribuzione_ogni_mese$Year[i], distribuzione_ogni_mese$Month[i], sep = "-")
    durata_evento <- 1
  } else if (distribuzione_ogni_mese$Colore_ssi[i] == "Rosso" && evento_in_corso) {
    durata_evento <- durata_evento + 1
  } else if (distribuzione_ogni_mese$Colore_ssi[i] == "Nero" && evento_in_corso) {
    evento_in_corso <- FALSE
    fine_evento <- paste(distribuzione_ogni_mese$Year[i-1], distribuzione_ogni_mese$Month[i-1], sep = "-")
    eventi_siccitosi <- rbind(eventi_siccitosi, data.frame(Data_di_Inizio = inizio_evento, Data_di_Fine = fine_evento, Durata = durata_evento))
    durata_evento <- 0
  }
}

# Se l'ultimo evento è ancora in corso alla fine del periodo, aggiungi l'evento
if (evento_in_corso) {
  fine_evento <- paste(distribuzione_ogni_mese$Year[nrow(distribuzione_ogni_mese)], distribuzione_ogni_mese$Month[nrow(distribuzione_ogni_mese)], sep = "-")
  eventi_siccitosi <- rbind(eventi_siccitosi, data.frame(Data_di_Inizio = inizio_evento, Data_di_Fine = fine_evento, Durata = durata_evento))
}

# Visualizza la nuova dataframe degli eventi siccitosi
print(eventi_siccitosi)


# Calcola la media delle durate
durata_media <- mean(eventi_siccitosi$Durata, na.rm = TRUE)

# Aggiungi la colonna "Durata Media" con la media delle durate
eventi_siccitosi$Durata_Media <- durata_media

# Visualizza il dataframe aggiornato
print(eventi_siccitosi)

# Inizializza una lista per memorizzare le somme di SSI per ogni evento
severity_list <- list()

# Calcola la somma di SSI per ciascun evento
for (i in 1:nrow(eventi_siccitosi)) {
  inizio_evento <- as.Date(paste0(eventi_siccitosi$Data_di_Inizio[i], "-01"))
  fine_evento <- as.Date(paste0(eventi_siccitosi$Data_di_Fine[i], "-01"))
  severity <- sum(distribuzione_ogni_mese$SSI[distribuzione_ogni_mese$Year * 100 + distribuzione_ogni_mese$Month >= as.numeric(format(inizio_evento, "%Y%m")) &
                                                distribuzione_ogni_mese$Year * 100 + distribuzione_ogni_mese$Month <= as.numeric(format(fine_evento, "%Y%m"))], na.rm = TRUE)
  severity_list[[i]] <- severity
}

# Aggiungi la colonna "Severity" al dataframe degli eventi siccitosi
eventi_siccitosi$Severity <- unlist(severity_list)

# Visualizza il dataframe aggiornato
print(eventi_siccitosi)


# Calcola la media dei valori nella colonna "Severity"
severity_media <- mean(eventi_siccitosi$Severity, na.rm = TRUE)

# Aggiungi la colonna "Severity_Media" con il valore della media dei valori di "Severity" per tutte le righe
eventi_siccitosi$Severity_Media <- severity_media

# Visualizza il dataframe aggiornato
print(eventi_siccitosi)


# Calcola la colonna "Intensità" come rapporto tra "Severity" e "Durata"
eventi_siccitosi$Intensità <- eventi_siccitosi$Severity / eventi_siccitosi$Durata

# Visualizza il dataframe aggiornato
print(eventi_siccitosi)


# Calcola la media della colonna "Intensità"
intensità_media <- mean(eventi_siccitosi$Intensità, na.rm = TRUE)

# Aggiungi la colonna "Intensità Media" con il valore medio calcolato
eventi_siccitosi$Intensità_Media <- intensità_media

# Visualizza il dataframe aggiornato
print(eventi_siccitosi)


# Calcola il numero totale di eventi nel dataframe
numero_eventi <- nrow(eventi_siccitosi)

# Aggiungi la colonna "Numero Eventi" con il numero totale di eventi
eventi_siccitosi$Numero_Eventi <- numero_eventi

# Visualizza il dataframe aggiornato
print(eventi_siccitosi)







