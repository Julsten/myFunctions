
n_trial_exclusions <- function(df_list, total=TRUE){
  if(!is.list(df_list) & !is.vector(df_list)) stop("df_list must be a list or vector")
  
  indices <- seq_along(df_list)
  output_table <- data.frame("Index"=indices)  
  if(!is.null(names(df_list)))output_table$Name <- names(df_list)
  
  output_table$n_trials <- sapply(df_list, nrow)
  output_table$n_exclusions <- NA
  output_table$n_exclusions[-1] <- output_table$n_trials[-length(output_table$n_trials)]-output_table$n_trials[-1]
  if(total){
    # total exclusions
    output_table[nrow(output_table)+1,]<-output_table[nrow(output_table),]
    output_table$Name[nrow(output_table)] <- "Total"
    output_table$n_exclusions[nrow(output_table)] <- nrow(df_list[[1]])-nrow(df_list[[length(df_list)]])
  }
  output_table$pct_exclusions <- round(output_table$n_exclusions/(output_table$n_exclusions+output_table$n_trials)*100,4)
  cnames <- c("Index","Data", "Trials (n)", "Excluded (n)", "Excluded (%)")
  if(is.null(names(df_list))){ cnames <- cnames[-2]}
  colnames(output_table) <- cnames
  
  return(output_table[-1])
  
}


sample_table <- function(sampledata) {
  #sample data needs variables participant, sex, hand and age!
  n <- length(unique(sampledata$participant))
  female <- prop.table(table(sampledata$sex))["f"]*100
  right <- prop.table(table(sampledata$hand))["r"]*100
  mean_age <- mean(sampledata$age)
  sd_age <- sd(sampledata$age)
  range_age <- range(sampledata$age)
  sample_table <- tibble(n, "female[%]" = female, "right handed[%]" =right, "mean age" = mean_age, "SD age"=sd_age, "age range" =paste(range_age[1], "-",range_age[2],sep=""))
  return(sample_table)
}



anova_eps_table <- function(x, digits=5, scipen=4, small=TRUE, pandoc=TRUE){
  suppressMessages({
    require(afex)
    require(tidyverse)
    require(dplyr)
    require(tibble)
    
    # check if sphericity is possible (any factor with > 2 levels)
    sphericity_check <- max(x$anova_table$`num Df`) > 1
    
    aovsum <- summary(x)
    
    # tibble for sphericity tests
    spht <- as.vector(aovsum$sphericity.tests)
    effect <- attr(aovsum$sphericity.test, "dimnames")[[1]]
    sphtcnames <- attr(aovsum$sphericity.test, "dimnames")[[2]]
    spht.table <- tibble::as_tibble(matrix(spht, nrow=length(effect)))
    colnames(spht.table) <- sphtcnames
    spht.table <- tibble::add_column(spht.table, effect, .before=T)
    
    
    # tibble for sphericity corrections
    effect <- attr(aovsum$pval.adjustments, "dimnames")[[1]]
    sph_corrections <- tibble::as_tibble(cbind(aovsum$pval.adjustments, rownames=NULL))
    sph_corrections <- tibble::add_column(sph_corrections, effect, .before=T)
    
    
    # tibble for anova 
    aovtb <- aovsum$univariate.tests
    effect <- (attr(aovtb, "dimnames"))[[1]]
    cnames <- (attr(aovtb, "dimnames"))[[2]]
    anova.table <- tibble::as_tibble(matrix(as.vector(aovtb), nrow=length(effect)))
    colnames(anova.table)<-cnames
    anova.table <- tibble::add_column(anova.table, effect, .before=T)
    
    if(pandoc){sig_sign = "*"} else{sig_sign="\\*"}
    # eta_squared
    ges <-c(NA, x$anova_table$ges)
    
    # join tibbles if spherictiy is possible
    if(sphericity_check) {
      temp_tib <- dplyr::full_join(anova.table, spht.table %>% dplyr::select(effect, `p-value`)) %>% 
        dplyr::rename(`Mauchly's p` = `p-value`)
      full_result_tib <- dplyr::full_join(temp_tib, sph_corrections %>% dplyr::select(effect, `GG eps`, `Pr(>F[GG])`))
      full_result_tib <- full_result_tib %>%  add_column(ges)
      
      full_result_tib <- full_result_tib %>% dplyr::mutate(sig= dplyr::case_when(
        is.na(`Mauchly's p`)&`Pr(>F)`<0.05 ~ sig_sign,
        `Mauchly's p` > 0.05&`Pr(>F)`<0.05 ~sig_sign,
        `Mauchly's p` < 0.05&`Pr(>F[GG])`<0.05 ~sig_sign,
        TRUE ~""
        
      ))
      if(small){
        full_result_tib <- full_result_tib %>% select(-`Sum Sq`, -`Error SS`) %>% 
          mutate(df=paste0(`num Df`, ", ", `den Df`)) %>% 
          mutate(GG = case_when(
            `Mauchly's p` < 0.05 ~ `GG eps`,
            TRUE ~ NA_real_
          )) %>% 
          mutate(`Pr(>F)`=case_when(
            `Mauchly's p` < 0.05 ~ `Pr(>F[GG])`,
            TRUE ~ `Pr(>F)`
          )) %>% 
          select(effect, df, `F value`, `Pr(>F)`, ges, GG, sig) %>% filter(effect!="(Intercept)")
        colnames(full_result_tib) <- c("effect", "df", "F", "p", "gen. eta^2", "GG", "sig")
        
        
      }
      
    }
    else{
      full_result_tib <- anova.table %>% dplyr::mutate(
        sig= dplyr::case_when(
          `Pr(>F)`<0.05 ~ sig_sign,
          TRUE ~ ""
        ))
      full_result_tib <- full_result_tib %>%  add_column(ges)
      if(small){
        full_result_tib <- full_result_tib %>% select(-`Sum Sq`, -`Error SS`) %>% 
          mutate( df = paste0(`num Df`, ", ", `den Df`)) %>% 
          select(effect, df, `F value`, `Pr(>F)`, ges, sig) %>% filter(effect!="(Intercept)")
        colnames(full_result_tib) <- c("effect", "df", "F", "p", "gen. eta^2","sig")
      }
    }
    # exclude sum of squares
    
    full_result_tib <- full_result_tib %>% dplyr::mutate_if(is.numeric, round, digits)
    rlang::local_options(scipen=scipen)
    return(full_result_tib)
  })
}


emm_contrasts_kable <- function(x, caption=NULL){ 
  
  require(kable)
  require(kableExtra)
  
  cons <- rbind(x$contrasts)
  
  # attributes for printing
  famSize <- attr(cons, "misc")$famSize
  adj <- attr(cons, "misc")$adjust
  
  foot_note_title <- paste("Adjustment Method for ", famSize, " Tests:")
  
  out <- footnote(kable_styling(kable(cons, caption=caption, )), general=adj, general_title=foot_note_title, footnote_as_chunk = T)
  return(out)
}


anova_kable <- function(x, caption=NULL) {
  suppressMessages({
    require(kableExtra)
    
    
    kable_styling(kable(anova_eps_table(x), caption=caption))
  })}





posthoc_table <- function(emm_obj){
  if(!"emm_list" %in% class(emm_obj)) stop("argument is not emm_list")
  output <- data.frame(rbind(emm_obj$contrasts, adjust="none"))
  output$adj.p.value <- data.frame(rbind(emm_obj$contrasts, adjust="bonf"))$p.value
  colnames(output) <- c("Contrast", "Difference", "SE","df", "t", "p", "p (adj.)")
  return(output)
} 



darken_plot <- function(plot){
  require(ggplot2)
  require(RColorBrewer)
  require(wesanderson)
  require(ggdark)
  
  if("patchwork" %in% class(plot)){
    plot <- plot & dark_theme_classic() & scale_fill_manual(values=wes_palette("Darjeeling1")) & theme( text=element_text(size=18), legend.position = "bottom") 
  }
  else{
    plot <- plot + dark_theme_classic() + scale_fill_manual(values=wes_palette("Darjeeling1")) + theme( text=element_text(size=20))
  }
  
  return(plot)
}



# helper function to check if a file is opened https://stackoverflow.com/questions/20038732/how-to-check-a-file-is-opened-or-closed-in-r
file_opened <- function(path) {
  suppressWarnings(
    "try-error" %in% class(
      try(file(path, 
               open = "w"), 
          silent = TRUE
      )
    )
  )
}


plots_to_ppt <- function(plot, doc, path='../output/plots/presentation_plots.pptx'){
  #doc needs to be specified with officer::read_pptx()
  if (file_opened(path)) stop("Error: PowerPoint Presentation needs to be closed!")
  if (!("rpptx" %in% class(doc))) stop("Error: doc needs to be of class rpptx (requires officer package)")
  if (!"ggplot" %in% class(plot)) stop("Error: plot needs to be of class ggplot")
  require(officer)
  require(rvg)
  doc <- add_slide(doc, 'Title and Content', 'Office Theme')
  vec_plt <- dml(ggobj=plot) # vector graphic
  doc <- ph_with(doc, vec_plt, ph_location_type(type="body")) 
  print(doc, target = path)
}


within_pd_barplot <- function(data, 
                              id="participant", 
                              condition="fixation_type", 
                              dv="estimation_error", 
                              errorwidth=0.05, 
                              barwidth=0.5,
                              barborder="black",
                              barbordersize=1,
                              errorcolor="grey", 
                              err_size=1,
                              hline=FALSE,
                              hline_size = 1){
  data_summary <- aggregate(data[dv], data[condition], FUN=mean)
  t.formula <- as.formula(paste(dv, "Group.2", sep=" ~ "))
  est_diff_t <- t.test(t.formula, aggregate(data[dv], 
                                            list(data[[id]], 
                                                 data[[condition]]), FUN=mean), paired=T)
  est_diff_conf <- max(est_diff_t$conf.int)- min(est_diff_t$conf.int)
  plot <- ggplot(data_summary, aes(x=.data[[condition]], y=.data[[dv]])) + 
    geom_bar(stat="identity", color=barborder, size=barbordersize, aes(fill=.data[[condition]]), width = barwidth) 
  if(hline) {
    plot <- plot + geom_hline(aes(yintercept=0), size=hline_size)
  }
    plot <- plot + geom_errorbar(aes(x=.data[[condition]], y=.data[[dv]], ymin=.data[[dv]]-est_diff_conf/2, 
                      ymax =.data[[dv]]+est_diff_conf/2),width=errorwidth, color=errorcolor, size=err_size)
  return(plot)
}











grouped_within_pd_barplot<-function(data, 
                                    id="participant", 
                                    group, 
                                    condition="fixation_type", 
                                    dv="estimation_error", 
                                    errorwidth=0.05,
                                    errsize=1, 
                                    barwidth=0.5,
                                    barbordersize=1, 
                                    barborder="black",
                                    bardodge=0.9,
                                    errorcolor="grey",
                                    hline=FALSE, 
                                    hline_size=1){
  require(tidyverse)
  out <- data.frame()
  data <- data %>% mutate(across(c(group, condition, id), factor))
  for(i in seq_along(levels(data[[group]]))){
  t_data <- aggregate(data[[dv]], list(data[[id]], data[[group]], data[[condition]]), FUN=mean)
  
  t_data <- t_data[t_data[["Group.2"]]==levels(t_data[["Group.2"]])[i],]
  
  exclude <- t_data %>% group_by(Group.1, Group.3, .drop=F) %>% tally() %>% filter(n==0) %>% .$Group.1
  if(length(exclude)>0) t_data <- t_data[t_data$Group.1!=exclude,]
  
  tt <- t.test(formula=x~Group.3, data=t_data, paired=T)
  conf_int <- max(tt$conf.int) - min(tt$conf.int)
  summary <- aggregate(t_data[["x"]], list(t_data[["Group.3"]]), FUN=mean)
  out <- rbind(out, cbind(group=levels(data[[group]])[i],summary, conf_int))
  }
  colnames(out) <- c(group, condition, "mean", "conf_int")
  plot <- ggplot(data=out, aes(x=.data[[group]], y=mean, fill=.data[[condition]])) +
    geom_bar(color=barborder, size=barbordersize, position=position_dodge(bardodge), stat="identity", width=barwidth)
  
  if(hline){
    plot <- plot + geom_hline(aes(yintercept=0), size=hline_size)
  }
  
  plot <- plot +geom_errorbar(aes(ymin=mean-conf_int/2, ymax=mean+conf_int/2),color=errorcolor, 
                  position=position_dodge(bardodge), width=errorwidth, size=errsize)
  return(plot)
}





# helper function for output printing
format_small_result <- function(x, digits = 3){
  if(digits==3){
    if(x < 0.001) x <- " < .001"
    else x <- paste0("= ",sub("^(-?)0.", "\\1.", sprintf("%.3f", x)))
  }
  if(digits==2){
    if(x < 0.001) x <- "< .01"
    else x <- paste0("= ",sub("^(-?)0.", "\\1.", sprintf("%.2f", x)))
  }
  return(x)
}

rm_anova_print <- function(aov_obj, effect) {
  require(knitr)
  sphericity_check <- max(aov_obj$anova_table$`num Df`) > 1
  aov_eps_tab <- anova_eps_table(aov_obj)
  
  F_val <- format(round(aov_eps_tab[aov_eps_tab$effect==effect,]$`F`,2), nsmall = 2)
  dfs <- aov_eps_tab[aov_eps_tab$effect==effect,]$df
  p_val <- format_small_result(aov_eps_tab[aov_eps_tab$effect==effect,]$p)
  eta <- format_small_result(round(aov_eps_tab[aov_eps_tab$effect==effect,]$`gen. eta^2` , 3))
  if(sphericity_check) {
    GG <- round(aov_eps_tab[aov_eps_tab$effect==effect,]$GG, 2)
    if(!is.na(GG)) GG <- format_small_result(GG, digits = 2)
  
    out <- paste0("*F*(", dfs, ")", " = ", F_val, ", *p* ", p_val, ", \u3b7~G~^2^ ", eta)
    if(!is.na(GG)) out <-  paste0(out, " (\u3b5 ", GG, ")")
  }
  else out <- paste0("*F*(", dfs, ")", " = ", F_val, ", *p* ", p_val, ", \u3b7~G~^2^ ", eta)
  knitr::asis_output(out)
}

# percentage format
pc_form <- function(x, digits=2) {
  x <- paste0(format(round(x, digits), nsmall=digits),"%")
  return(x)
}



#afex incomplete cases
afex_incomplete_cases <- function(x) {
  if(!("afex_aov" %in% class(x))) stop("x must by of class afex_aov!")
  return(attr(x$anova_table, "incomplete_cases") )
}


# fixing posthoc tables with leading X and no "ms" after miliseconds
fix_ms <- function(x) {
  
  return(paste(gsub("X", "", x),"ms"))
}



# read hdf5 files



# reading all condition variables from hdf5
read_condition_variables <- function(path="./data/hdf5/exp4_1.hdf5", 
                                     hdf5_path='data_collection/condition_variables/'){
  require(rhdf5)
  
  h5 <- H5Fopen(path)
  df_list <- h5&hdf5_path
  # list folder contents
  h5_ls <- h5ls(df_list)
  # valid condition_variable files (dim>0)
  val_dims <- h5_ls$dim != "0"
  names <- h5_ls$name[val_dims]
  
  df <- do.call(rbind, lapply(names, function(x) h5read(df_list, name=x)))
  
  return(df)
}




hdf5_to_data_frame <- function(path="./data/hdf5", hdf5_path='data_collection/condition_variables/'){
  path_list <- sapply(list.files(path, pattern="*.hdf5"), function(p) paste(path,p, sep="/"))
  df <- do.call(rbind, lapply(path_list, read_condition_variables,  hdf5_path = hdf5_path))
  return(df)
}



nice_plot <- function(gg_obj, 
                      title="Estimation Errors in Each Condition", 
                      caption=NULL, 
                      x_lab="Condition", 
                      y_lab=expression("Estimation Error [s] &  95%CI "[PD])) {
  out_plot <- gg_obj + scale_fill_manual(values=c("white", "black")) +  
    labs(x=x_lab, y=y_lab, title=title, caption=caption) +
    theme_classic()+
    theme(text=element_text(size=12), 
          axis.line = element_line(size=1), 
          axis.text = element_text(size=14),
          legend.position = "none", 
          axis.title.y = element_text( size=14, margin = margin(r = 0)),
          axis.title.x = element_text(margin = margin(t = 15)))
  return(out_plot)
}


# Test Data for grouped within paired difference barplot
# set.seed(123)
# a_1_mean <- 80
# a_2_mean <- 60
# b_1_mean <- 50
# b_2_mean <- 20
# a_1_sd <- 12
# a_2_sd <- 10
# b_1_sd <- 20
# b_2_sd <- 30
# 
# n_pp <- 40
# bias <- rnorm(n_pp, mean=10, sd=5)
# error <- function() rnorm(1, 0, 5)
# a_1 <- rnorm(n_pp, a_1_mean, a_1_sd)+bias
# a_2 <- rnorm(n_pp, a_2_mean, a_2_sd)+bias
# b_1 <- rnorm(n_pp, b_1_mean, b_1_sd)+bias
# b_2 <- rnorm(n_pp, b_2_mean, b_1_sd)+bias
# 
# df <- data.frame(vp=1:n_pp, a_1, a_2, b_1, b_2)
# library(tidyverse)
# df <- tidyr::pivot_longer(df, cols=!vp, names_to="cond") %>%
#   mutate(group=factor(case_when(cond %in% c("a_1", "a_2") ~"a", TRUE ~ "b"))) %>%
#   mutate(time = factor(case_when(cond %in% c("a_1", "b_1") ~ 1, TRUE ~ 2))) %>% select(-cond)


