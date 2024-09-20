
#about##########################################################################

  #Bayesian Hierarchical models (INLA) - MSOA level analyses for preterm birth

#setup##########################################################################


  #set location
  setwd("//ads.bris.ac.uk/Filestore/HealthSci SafeHaven/TOMMY/Ruta/code/Rmodels3")
  
  #clear environment/set print options
  rm(list=ls())
  options(max.print=99999)

  
#load libraries#################################################################

    
  library(rockchalk)
  library(rsample)
  library(dplyr)
  library(geojsonio)
  library(sf)
  library(spdep)
  library(spatialreg)
  library(INLA) 
  library(pROC)
  library(tmap)

    
#import data####################################################################

  
  #nmpa
  load("//ads.bris.ac.uk/filestore/HealthSci SafeHaven/TOMMY/Ruta/code/Rmodels3/dfprep.Rdata")
  #str(dfprep, list.len=ncol(dfprep))
  

#import geo boundaries and points###############################################
  
  
  #lad2017 boundaries
  sp_b_lad17<-geojson_read("Local_Authority_Districts_December_2017_Boundaries_GB_BUC.geojson",  what = "sp")
  summary(sp_b_lad17)
  sp_b_lad17 <- sp_b_lad17[substr(sp_b_lad17@data$LAD17CD,1,1) =="E", ] #leave england only
  
  #msoa2011 boundaries
  sp_b_msoa11<-geojson_read("MSOA_Dec_2011_Boundaries_Super_Generalised_Clipped_BSC_EW_V3.geojson",  what = "sp")
  summary(sp_b_msoa11)
  sp_b_msoa11 <- sp_b_msoa11[substr(sp_b_msoa11@data$MSOA11CD,1,1) =="E", ] #leave england only
  
  #msoa population weighted centroids
  sp_pwc_msoa11 <- geojson_read("MSOA_Dec_2011_PWC_in_England_and_Wales.geojson",  what = "sp")
  sp_pwc_msoa11 <- sp_pwc_msoa11[substr(sp_pwc_msoa11$msoa11cd,1,1) =="E", ] #leave england only
  
  #cities - boundaries
  sp_b_cities<-geojson_read("Major_Towns_and_Cities_Dec_2015_Boundaries_V2_2022.geojson",  what = "sp")
  summary(sp_b_cities)
  sp_b_cities <- sp_b_cities[sp_b_cities@data$TCITY15CD !="J01000064", ] #remove wales - newport
  sp_b_cities <- sp_b_cities[sp_b_cities@data$TCITY15CD !="J01000020", ] #remove wales - cardiff
  sp_b_cities <- sp_b_cities[sp_b_cities@data$TCITY15CD !="J01000098", ] #remove wales - swansea
  
  #regions - boundaries
  sp_b_rgn17<-geojson_read("Regions_December_2017_UGCB_in_England.geojson",  what = "sp")


#split sample testing###########################################################


  # Check if the column exists and remove it
  if ("split" %in% colnames(dfprep)) {
    dfprep <- dfprep[, !colnames(dfprep) %in% "split"]
  }

  set.seed(123)
  split <- initial_split(dfprep, prop = 0.80, strata = "outc") 
  dfprep1  <- training(split) #training sample 80%
  dfprep2  <- testing(split)  #internal validation sample 20%
  rm(split)

  #set split indicator
  dfprep1$split<-1
  dfprep2$split<-0

  #merge imputed into one - full dataset with indicator of split into test and validation
  dfprep <- dfprep1 %>%
  bind_rows(dfprep2)

  #clean up memory and files
  rm(dfprep1)
  rm(dfprep2)
  gc()


#set outcomes###################################################################

  
  #factor to numeric unmasked (for model fit - both validation and training sets included) 
    dfprep$outcn<-as.numeric(dfprep$outc)-1
    #table(dfprep$outcn,useNA= "always")  

  #factor to numeric masked (for training/validation - validation set is masked (hidden))
    dfprep$outcn_masked<-ifelse(dfprep$split==1,dfprep$outcn,NA)
    #table(dfprep$outcn_masked, useNA = "always")
  

#set formulas###################################################################


  #null model unmasked
    
    formula_0<-outcn~1
  
  #full model unmasked
  
    # Columns to exclude from the independent variables
    exclude_columns <- c("outc","outcn_masked", "LAD17CD", "MSOA11CD", "trustcode_cat","split")
    
    # Get the independent variables by excluding the specified columns
    independent_vars <- setdiff(colnames(dfprep), c("outcn", exclude_columns))
    
    # Create the formula object
    formula_f <- as.formula(paste("outcn ~", paste(independent_vars, collapse = " + ")))
    
    # Print the formula object
    print(formula_f)

  #null model masked
    
    formula_0m<-outcn_masked~ 1 

  
  #full model masked
    
    # Columns to exclude from the independent variables
    exclude_columns <- c("outc","outcn", "LAD17CD", "MSOA11CD", "trustcode_cat","split")
    
    # Get the independent variables by excluding the specified columns
    independent_vars <- setdiff(colnames(dfprep), c("outcn_masked", exclude_columns))
    
    # Create the formula object
    formula_fm <- as.formula(paste("outcn_masked ~", paste(independent_vars, collapse = " + ")))
    
    # Print the formula object
    print(formula_fm)
    
  #clean up
    
    rm(exclude_columns)
    rm(independent_vars)


#create neighbourhood matrices (MSOA)###########################################   


  #create sf object out of polygons and centroids data frames
  
    sp_b_msoa11_sf<-st_as_sf(sp_b_msoa11)
    #tm_shape(sp_b_msoa11_sf) +  tm_borders(lwd = 0.3, alpha = 0.4) 
    if (!all(st_is_valid(sp_b_msoa11_sf)))
      sp_b_msoa11_sf <- st_make_valid(sp_b_msoa11_sf)
    
    sp_pwc_msoa11_sf<-st_as_sf(sp_pwc_msoa11)

  #create adjacency based matrix (binary)
    
    sp_b_msoa11_sf |> poly2nb(queen = TRUE) -> sp_b_msoa11a
    sp_b_msoa11a
    
    matrix_adj_msoa <- spdep::nb2listw(sp_b_msoa11a, zero.policy = TRUE,style = "B")
    
    matrix_adj_msoa<- as(matrix_adj_msoa, "CsparseMatrix")
  
  #create knn inverse distance weighted matrix (symmetrical, row standardised)

    sp_pwc_msoa11_sf |> st_geometry() |>  st_centroid(of_largest_polygon = TRUE) -> centroids_msoa #centroids
    
    #distance2 - knn asymmetric
    ((centroids_msoa |> knearneigh(k=20) -> knn_k6) |> knn2nb() -> sp_b_msoa11c2)
    (sp_b_msoa11c2 |> n.comp.nb())$nc
    
    #distance3 - knn symmetric
    (knn_k6 |> knn2nb(sym = TRUE) -> sp_b_msoa11c3)
    (sp_b_msoa11c3 |> n.comp.nb())$nc


    sp_b_msoa11c3 |> 
      nbdists(centroids_msoa) |> 
      lapply(function(x) 1/(x/1000)) -> gwts #convert to kilometers and take inverse ???? 
    #lapply(function(x) 1/(x)) -> gwts #leave in meters and take inverse (too small!!!)
    (sp_b_msoa11c3 |> nb2listw(glist=gwts, style="W",zero.policy=TRUE) -> matrix_dist_msoa) |> 
      spweights.constants(zero.policy=TRUE) |> 
      data.frame() |> 
      subset(select=c(n, S0, S1, S2))

    matrix_dist_msoa<- as(matrix_dist_msoa, "CsparseMatrix")

  #clean up
    
    rm(sp_b_msoa11a)
    rm(sp_b_msoa11c3)
    rm(sp_b_msoa11c2)
    rm(centroids_msoa)
    rm(gwts)
    rm(knn_k6)

#set id for MSOA and masked/unmasked############################################

    ID2 <- as.integer(as.factor(dfprep$MSOA11CD))   
    ii <- which(is.na(dfprep$outcn_masked)==TRUE)
    i <- which(is.na(dfprep$outcn_masked)==FALSE)


#run INLA (MSOA) - null/full models with iidre/ssre (masked)###################


    #null model - iidre (unadjusted, multilevel, no spatial dependencies)

      gc()
      inla_iidre_msoa_0m <- inla(update(formula_0m, . ~ . + f(MSOA11CD, model = "iid")), family = "binomial", data = dfprep,
                                control.compute = list(dic = TRUE, waic = TRUE, cpo = TRUE))

    #null model - ssre (unadjusted, multilevel, with spatial dependencies)
      
      gc()
      inla_ssre_msoa_0m <- inla(update(formula_0m, . ~ . + f(ID2, model = "besag",graph = matrix_dist_msoa)), family = "binomial",  data = dfprep,
                               control.compute = list(dic = TRUE, waic = TRUE, cpo = TRUE,return.marginals.predictor = TRUE))


      
    #full model - iidre (adjusted, multilevel, no spatial dependencies)
      
      gc()
      inla_iidre_msoa_m <- inla(update(formula_fm, . ~ . + f(MSOA11CD, model = "iid")), family = "binomial", data = dfprep,
                        control.compute = list(dic = TRUE, waic = TRUE, cpo = TRUE))


    #full model - iidre (adjusted, multilevel, with spatial dependencies)


      gc()
      inla_ssre_msoa_m <- inla(update(formula_fm, . ~ . + f(ID2, model = "besag",graph = matrix_dist_msoa)), family = "binomial",  data = dfprep,
                        control.compute = list(dic = TRUE, waic = TRUE, cpo = TRUE,return.marginals.predictor = TRUE))



#run INLA (MSOA)  - null/full models with iidre/ssre (unmasked)#################
    
    
    #null model - iidre (unadjusted, multilevel, no spatial dependencies)
    
    gc()
    inla_iidre_msoa_0 <- inla(update(formula_0, . ~ . + f(MSOA11CD, model = "iid")), family = "binomial", data = dfprep,
                               control.compute = list(dic = TRUE, waic = TRUE, cpo = TRUE))
    
    #null model - ssre (unadjusted, multilevel, with spatial dependencies)
    
    gc()
    inla_ssre_msoa_0 <- inla(update(formula_0, . ~ . + f(ID2, model = "besag",graph = matrix_dist_msoa)), family = "binomial",  data = dfprep,
                              control.compute = list(dic = TRUE, waic = TRUE, cpo = TRUE,return.marginals.predictor = TRUE))
    
    
    
    #full model - iidre (adjusted, multilevel, no spatial dependencies)
    
    gc()
    inla_iidre_msoa <- inla(update(formula_f, . ~ . + f(MSOA11CD, model = "iid")), family = "binomial", data = dfprep,
                              control.compute = list(dic = TRUE, waic = TRUE, cpo = TRUE))
    
    
    #full model - iidre (adjusted, multilevel, with spatial dependencies)
    
    
    gc()
    inla_ssre_msoa <- inla(update(formula_f, . ~ . + f(ID2, model = "besag",graph = matrix_dist_msoa)), family = "binomial",  data = dfprep,
                             control.compute = list(dic = TRUE, waic = TRUE, cpo = TRUE,return.marginals.predictor = TRUE))
    

    
#save predictions in validation set############################################# 
    
    
    # List of models and their names
    models <- list(inla_iidre_msoa_0m = inla_iidre_msoa_0m, inla_iidre_msoa_m = inla_iidre_msoa_m, inla_ssre_msoa_0m = inla_ssre_msoa_0m, inla_ssre_msoa_m = inla_ssre_msoa_m)
    model_names <- names(models)
    
    # Initialize an empty data frame to store results
    results_pred <- data.frame(modelname_masked = character(), auc = numeric(), lowerci = numeric(), upperci = numeric(), 
                               tpr_at_10_fpr = numeric(), tpr_at_25_fpr = numeric(), stringsAsFactors = FALSE)
    
    
    # Original outcome values (validation set)
    original_values <- dfprep$outcn[ii]
    
    # Loop through each model to calculate AUC and 95% CI
    for (i in seq_along(models)) {
      # Get model and predicted values
      model <- models[[i]]
      predicted_values <- model$summary.fitted.values$mean[ii]
      
      # Calculate ROC, AUC, and CI
      roc_curve <- roc(original_values, predicted_values)
      auc_value <- auc(roc_curve)
      ci_value <- ci.auc(roc_curve)
      
      # Calculate TPR at 10% and 25% FPR
      tpr_at_10_fpr <- coords(roc_curve, x = 0.10, input = "fpr", ret = "tpr")
      tpr_at_25_fpr <- coords(roc_curve, x = 0.25, input = "fpr", ret = "tpr")
      
      
      # Store the results in the data frame
      results_pred <- rbind(results_pred, data.frame(
        modelname_masked = model_names[i],
        auc = auc_value,
        lowerci = ci_value[1],
        upperci = ci_value[3],
        tpr_at_10_fpr = as.numeric(tpr_at_10_fpr),
        tpr_at_25_fpr = as.numeric(tpr_at_25_fpr)
      ))
    }
    
    # Print the final results
    print(results_pred)
    
#save model fit and RE variance#################################################    
    
    # List of models and their names
    models <- list(inla_iidre_msoa_0 = inla_iidre_msoa_0, inla_iidre_msoa = inla_iidre_msoa, inla_ssre_msoa_0 = inla_ssre_msoa_0, inla_ssre_msoa = inla_ssre_msoa)
    model_names <- names(models)
    
    # Initialize an empty data frame to store results
    results_fit<- data.frame(modelname_unmasked = character(), cpo = numeric(), dic = numeric(), waic = numeric(), mlik = numeric(), 
                             var = numeric(), varupperci = numeric(), varlowerci = numeric(), stringsAsFactors = FALSE)
    

    # Loop through each model to calculate cpo,waic,mlik,dic, and random effect variance
    for (i in seq_along(models)) {
      # Get model and predicted values
      model <- models[[i]]

      cpo<- (-sum(log(model$cpo$cpo))) 
      dic<-model$dic$dic
      waic<-model$waic$waic
      mlik<-model$mlik
      var<- 1/model$summary.hyperpar$mean
      varupperci<- 1/model$summary.hyperpar$'0.025quant'
      varlowerci<- 1/model$summary.hyperpar$'0.975quant'
      
      # Store the results in the data frame
      results_fit <- rbind(results_fit, data.frame(
        modelname_unmasked = model_names[i],
        cpo = cpo,
        dic = dic,
        waic = waic,
        mlik = mlik[1],
        var= var,
        varlowerci= varlowerci,
        varupperci= varupperci

      ))
    }
    
    # Print the final results
    print(results_fit)
    
#merge results - predictive performance and fit#################################
    
    results_pred$row_id <- 1:nrow(results_pred)
    results_fit$row_id <- 1:nrow(results_fit)
    inla_results<- merge(results_pred, results_fit, by = "row_id")
    print(inla_results)

#export results#################################################################
    
    dir <- "//ads.bris.ac.uk/Filestore/HealthSci SafeHaven/TOMMY/Ruta/code/Rmodels3/output/"
    write.csv(inla_results, file = paste0(dir, "inla_results.csv"))
    
#save ssre msoa level results###################################################

    
    #update ID for plotting
      dfprep$ID2 <- as.integer(as.factor(dfprep$MSOA11CD))
      a<-subset(dfprep, select = c(MSOA11CD,ID2))
      a<-unique(a)
      a<-subset(a, select = c(MSOA11CD,ID2))
      if ("ID2" %in% colnames(sp_b_msoa11_sf)) { sp_b_msoa11_sf <- sp_b_msoa11_sf[, !colnames(sp_b_msoa11_sf) %in% "ID2"] }
      sp_b_msoa11_sf<-merge(sp_b_msoa11_sf,  a, by.x = "MSOA11CD", by.y = "MSOA11CD", all.x=TRUE)
      rm(a)
    
    #save random effect from null ssre model to sf object
      a<-inla_ssre_msoa_0$summary.random$ID2
      names(a)[names(a)=='0.025quant'] <- 'lower'
      names(a)[names(a)=='0.975quant'] <- 'upper'
      a$int <- (a$upper-a$lower)
      a$sign1<-ifelse(a$upper>0,ifelse(a$lower>0,1,0),0)
      a$sign2<-ifelse(a$upper<0,ifelse(a$lower<0,-1,0),0)
      a$sign<-ifelse(a$sign1==1,"Positive",ifelse(a$sign2==-1,"Negative",NA)) 
      a<-subset(a, select=c(ID,mean,int,sign))
      names(a)[names(a)=='mean'] <- 'inla0_ssre'
      names(a)[names(a)=='int'] <- 'inla0_ssre_int'
      names(a)[names(a)=='sign'] <- 'inla0_ssre_sign'
      if ("inla0_ssre" %in% colnames(sp_b_msoa11_sf)) { sp_b_msoa11_sf <- sp_b_msoa11_sf[, !colnames(sp_b_msoa11_sf) %in% "inla0_ssre"] }
      if ("inla0_ssre_int" %in% colnames(sp_b_msoa11_sf)) { sp_b_msoa11_sf <- sp_b_msoa11_sf[, !colnames(sp_b_msoa11_sf) %in% "inla0_ssre_int"] }
      if ("inla0_ssre_sign" %in% colnames(sp_b_msoa11_sf)) { sp_b_msoa11_sf <- sp_b_msoa11_sf[, !colnames(sp_b_msoa11_sf) %in% "inla0_ssre_sign"] }
      sp_b_msoa11_sf<-merge(sp_b_msoa11_sf,  a, by.x = "ID2", by.y = "ID", all.x=TRUE) 
      table(a$inla0_ssre_sign)
      rm(a)

#plot ssre msoa level results###################################################  

    #random effects - point estimates   
    plot_ssre<-tm_shape(sp_b_rgn17) +
      tm_borders(lwd = 0.3, alpha = 0.9)+
      tm_fill("grey", first=TRUE)+
      tm_shape(sp_b_msoa11_sf) +
      tm_fill(c("inla0_ssre"),midpoint = 0, colorNA = "grey", title = "SSRE", palette = "-RdBu",breaks = c(-Inf,-0.4,-0.3,-0.2,-0.1,0,0.1,0.2,0.3,0.4,+Inf))+ 
      tm_layout(legend.outside=F, legend.position = c("left","top"), panel.labels = c("England - preterm birth - unadjusted INLA (MSOA)"),
                panel.label.size=1.5, legend.text.size=1.3,legend.title.size=1.5) +
      tm_shape(sp_b_rgn17) +
      tm_borders(lwd = 0.3, alpha = 0.9) #+
    
    #random effects with p < 0.05
    plot_ssre_sign<-tm_shape(sp_b_rgn17) +
      tm_borders(lwd = 0.3, alpha = 0.9)+
      tm_fill("grey", first=TRUE)+
      tm_shape(sp_b_msoa11_sf) +
      tm_fill(c("inla0_ssre_sign"),midpoint = 0.00000000000000000001, colorNA = "grey90", textNA="Not \"interesting\"",  title = "SSRE (p≤0.05)", 
              palette = "-RdBu")+#,breaks = c(-1,1))+ #+ #breaks = c(-Inf,-0.25,-0.1,-0.05,0,0.05,0.1,0.25,+Inf)) +  
      tm_layout(legend.outside=F, legend.position = c("left","top"), panel.labels = c("England - preterm birth - undjusted INLA (MSOA)"),
                panel.label.size=1.5, legend.text.size=1.3,legend.title.size=1.5) +
      tm_shape(sp_b_rgn17) +
      tm_borders(lwd = 0.3, alpha = 0.9) 
    
    #random effect 95%CI width
    plot_ssre_int<-tm_shape(sp_b_rgn17) +
      tm_borders(lwd = 0.3, alpha = 0.9)+
      tm_fill("grey", first=TRUE)+
      tm_shape(sp_b_msoa11_sf) +
      tm_fill(c("inla0_ssre_int"),  title = "95%CI (SSRE)", 
              palette = "Blues",breaks = c(-Inf,0.45,0.5,0.55,0.6,0.65,0.7,+Inf)) +  
      tm_layout(legend.outside=F, legend.position = c("left","top"), panel.labels = c("England - preterm birth - unadjusted INLA (MSOA)"),
                panel.label.size=1.5, legend.text.size=1.3,legend.title.size=1.5) +
      tm_shape(sp_b_rgn17) +
      tm_borders(lwd = 0.3, alpha = 0.9) 

    
#export plots###################################################################
    

    exportplots<- function(file_name, plot_name)
    {png(file_name, res=100, height=730, width=650)
      print(plot_name)
      dev.off()}
    
    plot_ssre
    plot_ssre <- recordPlot()
    exportplots("//ads.bris.ac.uk/Filestore/HealthSci SafeHaven/TOMMY/Ruta/code/Rmodels3/output/plot_ssre_fig.png", plot_ssre)
    
    plot_ssre_sign
    plot_ssre_sign <- recordPlot()
    exportplots("//ads.bris.ac.uk/Filestore/HealthSci SafeHaven/TOMMY/Ruta/code/Rmodels3/output/plot_ssre_sign_fig.png", plot_ssre_sign)
    
    plot_ssre_int
    plot_ssre_int <- recordPlot()
    exportplots("//ads.bris.ac.uk/Filestore/HealthSci SafeHaven/TOMMY/Ruta/code/Rmodels3/output/plot_ssre_int_fig.png", plot_ssre_int)
    
    
    
#save workspace#############################################################
    
    rm(list=setdiff(ls(), c("dfprep", 
                            "inla_iidre_msoa","inla_iidre_msoa_0","inla_iidre_msoa_0m","inla_iidre_msoa_m",
                            "inla_results","inla_ssre_msoa","inla_ssre_msoa_0","inla_ssre_msoa_0m","inla_ssre_msoa_m"
    )))
    
    
    save.image(file = "//ads.bris.ac.uk/Filestore/HealthSci SafeHaven/TOMMY/Ruta/code/Rmodels3/inla_models.RData")
    
