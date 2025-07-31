# Load required libraries
suppressMessages(library("dplyr"))
suppressMessages(library("lubridate"))
suppressMessages(library("sf"))
# suppressMessages(library("ggplot2"))

# Get the script path and set working directory
script_path <- commandArgs(trailingOnly=FALSE)
script_path <- script_path[grep("--file=", script_path)]
script_path <- sub("--file=", "", script_path)
script_dir <- ifelse(length(script_path) > 0,
                     dirname(script_path),
                     dirname(rstudioapi::getSourceEditorContext()$path))
setwd(dirname(script_dir))

# LOG - Print the current working directory
cat("R Current directory:", getwd(), "\n")

# Import raw shapefiles
cat("Importing shapefiles... ")
reg_sf <- read_sf("raw/Limiti01012020/Reg01012020/Reg01012020_WGS84.shp")
prov_sf <- read_sf("raw/Limiti01012020/ProvCM01012020/ProvCM01012020_WGS84.shp")
mun_sf <- read_sf("raw/Limiti01012020/Com01012020/Com01012020_WGS84.shp")
cat("Done!\n")

# Import raw data from all sources
cat("Importing raw data... ")
mortality_df_raw <- read.csv("raw/comuni_giornaliero_31dicembre.csv", encoding = "latin1")
population_df_raw <- read.csv("raw/POSAS_2021_it_Comuni.csv", sep = ";", skip = 1)
income_df_names <- names(read.csv("raw/Redditi_e_principali_variabili_IRPEF_su_base_comunale_CSV_2020.csv", sep=";", header = T, nrows = 0))
income_df_raw <- read.csv("raw/Redditi_e_principali_variabili_IRPEF_su_base_comunale_CSV_2020.csv", sep = ";", skip = 1, header = F) %>%
  select(-c("V51","V52")) %>% `colnames<-`(income_df_names) %>% filter(Regione != "Mancante/errata")
rm(income_df_names)
cat("Done!\n")

# Merge Montericcardo into Pesaro municipality
cat("Merging Monteciccardo into Pesaro municipality... ")
merged_mun <- mun_sf %>%
  filter(COMUNE %in% c("Monteciccardo", "Pesaro")) %>%
  summarise(COD_RIP = 3, COD_REG = 11, COD_PROV = 41, COD_CM = 0, COD_UTS = 41,
    PRO_COM = 41044, PRO_COM_T = "041044", COMUNE = "Pesaro", COMUNE_A = NA, CC_UTS = 1,
    SHAPE_LENG = as.numeric(st_length(st_boundary(st_union(geometry)))),
    SHAPE_AREA = as.numeric(st_area(st_union(geometry))),
    geometry = st_union(geometry)
  )
mun_sf <- mun_sf %>% filter(!COMUNE %in% c("Monteciccardo", "Pesaro")) %>% bind_rows(merged_mun)
rm(merged_mun)
cat("Done!\n")

# Data clean-up to extract useful information
cat("Cleaning and processing data... ")
mortality_df <- mortality_df_raw %>%
  mutate(ME = as.factor(month(as.Date(sprintf("2024%04g", GE), format="%Y%m%d")))) %>%
  rename(COD_REG = REG, COD_PROV = PROV) %>%  filter(ME %in% c("3","4","5")) %>%
  select(-c("TIPO_COMUNE", "GE", matches("(M|F)_[0-9+]"))) %>%
  group_by(COD_REG, COD_PROV, NOME_REGIONE, NOME_PROVINCIA, NOME_COMUNE, COD_PROVCOM) %>%
  summarise(across(matches("T_[0-9+]"), sum), .groups = "drop") %>%
  rowwise() %>% mutate(ET = mean(c_across(T_11:T_19), na.rm = TRUE)) %>% ungroup() %>%
  select(!matches("T_1[1-9]"))
cat("Done!\n")

# Computing the % ofover-65 population and population density in each municipality
cat("Computing population statistics... ")
population_df <- population_df_raw %>%
  filter(Età != 999) %>%
  group_by(Codice.comune, Comune) %>%
  summarise(TotalPop = sum(Totale), Over65 = sum(if_else(Età >= 65, Totale, 0)), .groups = "drop") %>%
  mutate(PercOver65 = Over65 / TotalPop) %>%
  left_join(mun_sf, by = c("Codice.comune" = "PRO_COM")) %>%
  mutate(PopDens = TotalPop / (SHAPE_AREA * 1e-6)) %>%
  select(Codice.comune, Comune, TotalPop, PercOver65, PopDens)
cat("Done!\n")

# Compute Average Personal Income in each municipality
cat("Computing income statistics... ")
income_df <- income_df_raw %>%
  select(!matches("Reddito.[^c]|Imposta|Bonus|Addizionale|Reddito.complessivo.m")) %>% rowwise() %>%
  mutate(TotFreq = sum(c_across(matches("Frequenza")), na.rm = TRUE),
         TotIncome = sum(c_across(matches("Ammontare.in.euro")), na.rm = TRUE),
         MeanIncome = TotIncome / TotFreq) %>% ungroup() %>%
  select(Codice.Istat.Comune, Denominazione.Comune, TotFreq, TotIncome, MeanIncome)
cat("Done!\n")

# Merging data from all sources
cat("Merging all datasets... ")
real_data_sf <- mun_sf %>%
  left_join(mortality_df, by = c("PRO_COM" = "COD_PROVCOM", "COD_REG", "COD_PROV")) %>%
  left_join(population_df, by = c("PRO_COM" = "Codice.comune")) %>%
  left_join(income_df, by = c("PRO_COM" = "Codice.Istat.Comune")) %>%
  select(COD_REG, NOME_REGIONE, COD_PROV, NOME_PROVINCIA, PRO_COM, COMUNE, T_20, ET, TotalPop, PercOver65, PopDens, TotFreq, TotIncome, MeanIncome) %>%
  rename(COD_MUN = PRO_COM, NAME_REG = NOME_REGIONE, NAME_PROV = NOME_PROVINCIA, NAME_MUN = COMUNE) %>%
  mutate(TotIncome = as.character(TotIncome))
cat("Done!\n")

# Fill-up info from two municipalities with missing data
to_fill <- which(is.na(real_data_sf$ET))
real_data_sf[to_fill, "NAME_REG"] <- "Piemonte"
real_data_sf[to_fill, "NAME_PROV"] <- c("Cuneo","Vercelli")

# Save final dataset for analysis
dir.create(file.path("input"), recursive = T)
st_write(real_data_sf, file.path("input","full_dataset.shp"), delete_layer = TRUE)
