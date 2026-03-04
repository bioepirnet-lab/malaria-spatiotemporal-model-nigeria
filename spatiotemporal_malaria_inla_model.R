library(INLA)
library(sf)
library(spdep)
library(tidyverse)
library(haven)
library(stringr)
library(dplyr)
library(tidyr)


gamdata <- read_sav("MergeDataEdit.SAV")
gamdata=as.data.frame(gamdata)
gamdata=gamdata[,c("HHID","HV025","HV040","bednetUse","HV270","HC1","HC27","HC56","HC57",
                   "HC61","HML32","yearstudy","agecat","SHSTATE")]## select variables from gamdata 
view(gamdata)

gamdata=na.omit(gamdata)
nrow(gamdata)

################## Data preparation ############

gamdata$HV025=factor(gamdata$HV025,levels = c(1,2),labels= c("Urban","Rural"))
gamdata$bednetUse=factor(gamdata$bednetUse,levels = c(0,1),
                         labels =c("No","Yes") )
gamdata$HV270=factor(gamdata$HV270,levels = c(1,2,3,4,5),
                     labels =c("Poorest","Poorer","Middle","Richer",
                               "Richest") )
gamdata$HC57=factor(gamdata$HC57,levels = c(1,2,3,4),
                    labels =c("Severe","Moderate","Mild","NonAnemic") )
gamdata$HC61=factor(gamdata$HC61,levels = c(0,1,2,3),
                    labels =c("noEdu","Primary","Secondary","Higher") )
gamdata$HC27=factor(gamdata$HC27,levels = c(1,2),labels = c("Male","Female"))
gamdata <- filter(gamdata, HML32 %in% c(0,1))
gamdata$yearstudy=factor(gamdata$yearstudy,levels = c(1,2,3),
                         labels = c("2010","2015","2021"))

gamdata$agecat=factor(gamdata$agecat,levels = c(1,2,3),
                      labels = c("Infants","Toddlers","Preschoolars"))
gamdata$SHSTATE <- gamdata$SHSTATE

view(gamdata)
reg_data <- data.frame(
  REGCODE = c(10, 20, 30, 40, 50, 60, 70, 80, 90, 100, 110, 120, 130, 140, 150, 160, 170, 180, 190, 200, 210, 220, 230, 240, 250, 260, 270, 280, 290, 300, 310, 320, 330, 340, 350, 360, 370),
  REGNAME = c('Sokoto', 'Zamfara', 'Katsina', 'Jigawa', 'Yobe', 'Borno', 'Adamawa', 'Gombe', 'Bauchi', 'Kano', 'Kaduna', 'Kebbi', 'Niger', 'Federal Capital Territory', 'Nasarawa', 'Plateau', 'Taraba', 'Benue', 'Kogi', 'Kwara', 'Oyo', 'Osun', 'Ekiti', 'Ondo', 'Edo', 'Anambra', 'Enugu', 'Ebonyi', 'Cross River', 'Akwa Ibom', 'Abia', 'Imo', 'Rivers', 'Bayelsa', 'Delta', 'Lagos', 'Ogun')
)

gamdata$SHSTATE <- as.numeric(gamdata$SHSTATE)
reg_data$REGCODE <- as.numeric(reg_data$REGCODE)
gamdata <- left_join(gamdata, reg_data, by = c("SHSTATE" = "REGCODE"))
cNames = c("CASEID","residence","clusterAltitude","bednetUse","wealthIndex",
           "age","sex","hemoglobin","anemiaLevel","motherEdlevel",
           "microscopy","yearstudy","agecat","state", "REGNAME")
colnames(gamdata)=cNames

view(gamdata)


######### Read in the shape file

map <- st_read("C:/Users/DELL/Documents/SPATIAL MODELLING/NIGERIA MAP/gadm41_NGA_shp", layer = "gadm41_NGA_1")

### Create Adjacent matrix and the time index

nb<-poly2nb(map)
nb2INLA("NGR.graph", nb)
NGR.adj <- inla.read.graph(filename="NGR.graph")


dw <- reshape(gamdata,
              timevar = 'yearstudy',
              idvar = 'REGNAME',
              direction = "wide"
)


map <- merge(map, dw, by.x = "NAME_1", by.y = "REGNAME") 
map_sf <- st_as_sf(map)

if (is.factor(gamdata$yearstudy)) {
  gamdata$yearstudy <- as.numeric(as.character(gamdata$yearstudy))
}

gamdata$idtime <- 1 + gamdata$yearstudy - min(gamdata$yearstudy)
gamdata$state <- as.numeric(as.factor(gamdata$state))
gamdata$state_spatial <- gamdata$state 
gamdata$state_iid <- as.numeric(gamdata$state)


#### Write the INLA formula and the INLA call

formula <- microscopy ~residence + bednetUse + wealthIndex + sex + anemiaLevel+motherEdlevel+yearstudy+
  f(idtime, model = "rw1") +
  f(age, model = "rw1") + 
  f(state_spatial, model = "besag", graph = NGR.adj) + 
  f(state_iid, idtime, model = "iid")

mod <- inla(formula, data = gamdata, family = "binomial", 
            control.compute = list(dic = TRUE, waic = TRUE), verbose = TRUE)
summary(mod)

gamdata$R <- mod$summary.fitted.values[, "mean"]

View(gamdata)

### ### Aggregate fitted probabilities by state #########

state_year_pred <- gamdata %>%
  dplyr::group_by(REGNAME, yearstudy) %>%
  dplyr::summarise(
    R = mean(R, na.rm = TRUE),
    .groups = "drop"
  )

names(state_year_pred)

map_sf_states <- map_sf %>%
  dplyr::select(NAME_1, geometry) %>%
  distinct()

map_sf_plot <- map_sf_states %>%
  left_join(state_year_pred, by = c("NAME_1" = "REGNAME"))

###### Plot 

ggplot(map_sf_plot) +
  geom_sf(aes(fill = R)) +
  facet_wrap(~yearstudy) +
  scale_fill_distiller(
    palette = "RdYlBu",
    direction = -1,
    name = "Prevalence"
  ) +
  theme_minimal() +
  theme(
    legend.position = "right",
    strip.text = element_text(face = "bold")
  )


###################################################################
############  Bivariate Analysis######################### 
##################################################################


##### proportion Percentage  across fixed effect ######

fixed_vars <- c(
  "yearstudy",
  "residence",
  "wealthIndex",
  "motherEdlevel",
  "bednetUse",
  "sex",
  "anemiaLevel"
)


create_table_manuscript_style <- function(data, vars){
  
  bind_rows(
    lapply(vars, function(v){
      
      # Age breakdown within each category
      age_part <- data %>%
        group_by(.data[[v]], agecat) %>%
        summarise(
          N = n(),
          Pos = sum(microscopy == 1, na.rm = TRUE),
          .groups = "drop"
        )
      
      # Row totals (denominator for age %)
      totals <- age_part %>%
        group_by(.data[[v]]) %>%
        summarise(Row_total = sum(N), .groups = "drop")
      
      age_part <- age_part %>%
        left_join(totals, by = v) %>%
        mutate(
          N_percent = round(100 * N / Row_total),
          Pos_percent = round(100 * Pos / N),
          N_display = paste0(N, " (", N_percent, ")"),
          Pos_display = paste0(Pos, " (", Pos_percent, ")")
        ) %>%
        select(.data[[v]], agecat, N_display, Pos_display)
      
      # Pivot to wide
      age_wide <- age_part %>%
        pivot_wider(
          names_from = agecat,
          values_from = c(N_display, Pos_display)
        )
      
      # Total column
      total_part <- data %>%
        group_by(.data[[v]]) %>%
        summarise(
          N_total = n(),
          Pos_total = sum(microscopy == 1, na.rm = TRUE),
          .groups = "drop"
        ) %>%
        mutate(
          Pos_percent_total = round(100 * Pos_total / N_total),
          N_display_total = paste0(N_total, " (", 
                                   round(100 * N_total / nrow(data)), ")"),
          Pos_display_total = paste0(Pos_total, " (", Pos_percent_total, ")")
        )
      
      final_block <- age_wide %>%
        left_join(total_part, by = v) %>%
        rename(Category = !!v) %>%
        mutate(
          Category = as.character(Category),
          Variable = v,
          .before = 1
        )
      
      final_block
    })
  )
}

final_table <- create_table_manuscript_style(gamdata, fixed_vars)


View(final_table)

library(openxlsx)

write.xlsx(
  final_table,
  file = "Table1_Study_Characteristics.xlsx",
  sheetName = "Bivariate Table",
  rowNames = FALSE
)


####### proportion Percentage  across state #############

create_state_table <- function(data){
  
  overall_total <- nrow(data)
  
  # Age breakdown per state
  age_part <- data %>%
    group_by(REGNAME, agecat) %>%
    summarise(
      N = n(),
      Pos = sum(microscopy == 1, na.rm = TRUE),
      .groups = "drop"
    )
  
  # State totals (denominator for age %)
  state_totals <- age_part %>%
    group_by(REGNAME) %>%
    summarise(State_total = sum(N), .groups = "drop")
  
  age_part <- age_part %>%
    left_join(state_totals, by = "REGNAME") %>%
    mutate(
      N_percent = round(100 * N / State_total),
      Pos_percent = round(100 * Pos / N),
      N_display = paste0(N, " (", N_percent, ")"),
      Pos_display = paste0(Pos, " (", Pos_percent, ")")
    ) %>%
    select(REGNAME, agecat, N_display, Pos_display)
  
  # Pivot to wide format
  age_wide <- age_part %>%
    pivot_wider(
      names_from = agecat,
      values_from = c(N_display, Pos_display)
    )
  
  # Total columns
  total_part <- data %>%
    group_by(REGNAME) %>%
    summarise(
      N_total = n(),
      Pos_total = sum(microscopy == 1, na.rm = TRUE),
      .groups = "drop"
    ) %>%
    mutate(
      N_percent_total = round(100 * N_total / overall_total),
      Pos_percent_total = round(100 * Pos_total / N_total),
      N_display_total = paste0(N_total, " (", N_percent_total, ")"),
      Pos_display_total = paste0(Pos_total, " (", Pos_percent_total, ")")
    )
  
  final_state_table <- age_wide %>%
    left_join(total_part, by = "REGNAME") %>%
    rename(State = REGNAME) %>%
    arrange(State)
  
  return(final_state_table)
}


state_table <- create_state_table(gamdata)


write.xlsx(
  state_table,
  "Supplementary_Table_State_Level.xlsx",
  sheetName = "State Distribution",
  rowNames = FALSE
)

###############################################################

# --- Sensitivity Analysis across states ---

############################################################

orig_model <- inla(formula, family = "binomial", data = gamdata,
                   control.predictor = list(compute = TRUE))

# 1. Compute the original mean risk from the full model
orig_mean <- mean(mod$summary.fitted.values$mean, na.rm = TRUE)

# 2. Create an empty data frame
sensitivity_results <- data.frame(State = character(),
                                  OrigMean = numeric(),
                                  ExcludedMean = numeric(),
                                  Sensitivity = numeric(),
                                  stringsAsFactors = FALSE)

# 3. Loop through states
for (st in unique(gamdata$REGNAME)) {
  cat("Running LOO for:", st, "\n")
  
  excluded_data <- gamdata[gamdata$REGNAME != st, ]
  
  loo_mod <- inla(formula, family = "binomial", data = excluded_data,
                  control.predictor = list(compute = TRUE))
  
  excluded_mean <- mean(loo_mod$summary.fitted.values$mean, na.rm = TRUE)
  
  sens_val <- abs(excluded_mean - orig_mean) / orig_mean
  
  sensitivity_results <- rbind(sensitivity_results,
                               data.frame(State = st,
                                          OrigMean = orig_mean,
                                          ExcludedMean = excluded_mean,
                                          Sensitivity = sens_val))
}

# 4. Sort by sensitivity
sensitivity_results <- sensitivity_results[order(-sensitivity_results$Sensitivity), ]

# 5. View results
print(head(sensitivity_results, 10))

# 6. Export to CSV
write.csv(sensitivity_results, "sensitivity_results.csv", row.names = FALSE)

print(sensitivity_results)


##################################################################
############  Sensitivity Analysis for prior ###############
###################################################################

#### Set up: the baseline formula

build_formula <- function(pc_u = 1, pc_alpha = 0.01) {
  
  f_spatial <- sprintf(
    "f(state_spatial, model='besag', graph=NGR.adj,
        hyper=list(prec=list(prior='pc.prec', param=c(%f, %f))))",
    pc_u, pc_alpha
  )
  
  as.formula(paste0(
    "microscopy ~ residence + bednetUse + wealthIndex + sex + anemiaLevel +
     motherEdlevel + yearstudy +
     f(idtime, model='rw1') +
     f(age, model='rw1') + ",
    f_spatial, " +
     f(state_iid, idtime, model='iid')"
  ))
}

##### Function to run models with modifiable priors
run_inla_model <- function(pc_u = 1, fixed_prec = 1/1e6, pc_alpha = 0.01) {
  
  formula_mod <- build_formula(pc_u = pc_u, pc_alpha = pc_alpha)
  
  mod <- inla(
    formula_mod,
    data = gamdata,
    family = "binomial",
    control.fixed = list(mean = 0, prec = fixed_prec),
    control.compute = list(dic = TRUE, waic = TRUE),
    verbose = TRUE
  )
  
  return(mod)
}

###### Define the prior scenarios

sensitivity_scenarios <- list(
  baseline = list(fixed_prec = 1/1e6, pc_u = 1),
  fixed_mod = list(fixed_prec = 1/1e3, pc_u = 1),
  pc_stronger = list(fixed_prec = 1/1e6, pc_u = 0.2),
  pc_weaker = list(fixed_prec = 1/1e6, pc_u = 2)
)

#### Run all models automatically

results_list <- list()

for (name in names(sensitivity_scenarios)) {
  params <- sensitivity_scenarios[[name]]
  message("Running model: ", name)
  
  results_list[[name]] <- run_inla_model(
    pc_u = params$pc_u,
    fixed_prec = params$fixed_prec,
    pc_alpha = 0.01
  )
}

#### Extract summaries

extract_summary <- function(model, label) {
  
  # --- FIXED EFFECTS ---
  if (!is.null(model$summary.fixed)) {
    fixed <- cbind(
      Model = label,
      Effect = rownames(model$summary.fixed),
      as.data.frame(model$summary.fixed)
    )
  } else {
    fixed <- data.frame(Model=label, Effect=NA)
  }
  
  # --- HYPERPARAMETERS ---
  if (!is.null(model$summary.hyperpar)) {
    hyper <- cbind(
      Model = label,
      Parameter = rownames(model$summary.hyperpar),
      as.data.frame(model$summary.hyperpar)
    )
  } else {
    hyper <- data.frame(Model=label, Parameter=NA)
  }
  
  # --- FIT STATISTICS ---
  fit <- data.frame(
    Model = label,
    DIC = ifelse(!is.null(model$dic$dic), model$dic$dic, NA),
    WAIC = ifelse(!is.null(model$waic$waic), model$waic$waic, NA)
  )
  
  return(list(fixed = fixed, hyper = hyper, fit = fit))
}


### Apply extraction to all models

summaries <- lapply(names(results_list), function(n) {
  extract_summary(results_list[[n]], n)
})
names(summaries) <- names(results_list)
names(summaries) <- names(results_list)


#### Combine all results

fixed_effects_all <- do.call(rbind,
                             lapply(names(summaries), function(n) summaries[[n]]$fixed))

random_effects_all <- do.call(rbind,
                              lapply(names(summaries), function(n) summaries[[n]]$hyper))

model_fit_stats <- do.call(rbind,
                           lapply(names(summaries), function(n) summaries[[n]]$fit))

head(fixed_effects_all)
head(random_effects_all)
head(model_fit_stats)

print(fixed_effects_all)
print(random_effects_all)
print(model_fit_stats)




