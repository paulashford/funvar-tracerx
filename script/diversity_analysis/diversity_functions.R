# diversity_functions.R
# Data wrangling to convert FIE count exports to suitable iNEXT analysis formats
# for examples of calls see:
#   copy_number/analysis/07_diversity_iNEXT/07_diversity_simplified.R 
# 24 02 2023

# Initial filter of table rows eg by cancer type, FIE/muts etc...
filter_analysis <- function( df_analysis, 
                                type_dat = 'FIE', # 'nf', 'nd', 'nfd' 
                                variant_class = c( 'Missense_Mutation', 'Silent' ), 
                                txp = c( 'Y', 'N' ),
                                timing_in = c( 'early', 'late' ),
                                cancer_in = c( 'LUAD', 'LUSC' ),  # Or "NSCLC"
                                remove_cols = c( 'VARIANT_CLASS', 'SCORE_THRESHOLD', 'DRIVERMUT', 'MUTANT_EXPRESSED', 'NUM_ROWS' ),
                                simpler_labels = simplify_lables,
                                pyclone_filt = c( 'none' )  # eg for clonal only 
                            ){
    # Filter
    df_analysis <- df_analysis %>%
        ungroup() %>%
        filter(VARIANT_CLASS %in% variant_class ) %>% 
        filter( TXP %in% txp ) %>% 
        filter( OVERALL_TIMING %in% timing_in )
    
    # Clonality
    if ( pyclone_filt[[1]] != 'none' ){
        df_analysis <- df_analysis %>%
            ungroup() %>%
            filter( PYCLONECLONAL_SC %in% pyclone_filt )    
    }

    # Only filter cancer type if not "NSCLC"
    ctypes <- paste0( cancer_in, collapse = '-' )
    if ( ctypes != 'NSCLC' ){
            df_analysis <- filter( df_analysis, CANCER_TYPE %in% cancer_in )
        }
        
    # FIE or background?
        if ( type_dat == 'FIE' ){
            df_analysis <- filter( df_analysis, SCORE_THRESHOLD == 'MC3' ) 
        } else if ( type_dat == 'nf' ){
            df_analysis <- filter( df_analysis, is.na( SCORE_THRESHOLD ) )
        } else if ( type_dat == 'nd' ){
            df_analysis <- filter( df_analysis, is.na( DRIVERMUT ) )
        } else if ( type_dat == 'nfd' ){
            df_analysis <- filter( df_analysis, is.na( SCORE_THRESHOLD ) & is.na( DRIVERMUT ) )  
        } else{
            print( "invalid type_dat, use: FIE, 'nf', 'nd', 'nfd'" )
            return( -1 )
        }
        
    # Create combo column for TxP if both Y and N
    ptxp <- paste0( txp, collapse = '-' )
    if ( ( ptxp == 'N-Y' ) | ( ptxp == 'Y-N' ) ){
        ptxp <- 'Tx-TCGA'
        df_analysis <- mutate( df_analysis, TXP = ptxp )
    }
    
    # Create combo for cancer types
    # LUAD LUSC combo
    if ( ( ctypes == 'LUAD-LUSC' ) | ( ctypes == 'LUSC-LUAD' ) ){
        ctypes <- 'LUAD-LUSC'
        df_analysis <- mutate( df_analysis, CANCER_TYPE = ctypes )
    # this would include other NSCLC not just LUAD LUSC...
    }else if ( ctypes == 'NSCLC'){
        df_analysis <- mutate( df_analysis, CANCER_TYPE = ctypes )
    }

    # Add data type (FIE/nfd etc) to SPECIES_TYPE for plotting etc
    # First check if including Silent mutations
    varclasses <- paste0( variant_class, collapse = '-' ) 
    if ( ( varclasses == 'Missense_Mutation-Silent' ) | ( varclasses == 'Silent-Missense_Mutation' ) ){
        varclasses <- 'Missense_Mutation-Silent'
        df_analysis <- mutate( df_analysis, VARIANT_CLASS = varclasses )
    }
    # now update SPECIES TYPE to include variant classes (not if simpler_labels==TRUE)
    if ( simpler_labels == FALSE ){
        if ( varclasses == 'Missense_Mutation-Silent' ){
            df_analysis <- mutate( df_analysis, SPECIES_TYPE = paste0( 'mis-syn-', SPECIES_TYPE ) )
        }else if ( varclasses == 'Missense_Mutation' ) {
            df_analysis <- mutate( df_analysis, SPECIES_TYPE = paste0( 'mis-', SPECIES_TYPE ) )
        }else if ( varclasses == 'Silent' ) {
            df_analysis <- mutate( df_analysis, SPECIES_TYPE = paste0( 'syn-', SPECIES_TYPE ) )
        }
    }

    # Add in type_dat (FIE, nfd, or the subsets for FIEs and mutations together)
    # in future may update SPECIES_TYPE in groups (not implemented)
    if ( type_dat %in% c( 'FIE', 'nf', 'nd', 'nfd' ) ){
        # in general, just add this to SPECIES_TYPE if just FIES or just background
        # NB simpler labels for easier to read plots keys!
        type_lab <- type_dat
        if ( simpler_labels == TRUE ){
            if ( type_lab == 'nfd' ){ type_lab <- 'background mutations' }
            if ( type_lab == 'FIE' ){ type_lab <- 'FIEs' }
            df_analysis <- mutate( df_analysis, SPECIES_TYPE = paste0( type_lab ) )
        }else{
            df_analysis <- mutate( df_analysis, SPECIES_TYPE = paste0( type_lab, '-', SPECIES_TYPE ) )
        }
        
    }

    # Remove cols not req in select
    df_analysis <- df_analysis %>%
        select( -all_of(remove_cols) )
    
    # chuck back filtered dataset
    return( df_analysis )
}

# Groups filtered dataset  by selecting correct columns
# note that the 'gene', 'funfam', 'gene-funfam' makes no difference as these group
# under the species column
# creat tot_muts as sum(...) - necessary as the standard df_analysis format is broken down
# by gene, funfam, txp, variant_class etc and so will need to group/sum/distinct ready for iNEXT
group_analysis <- function( df_analysis, 
                            analysis_type = 'timing' # 'duplicated-clonality' 
                        ){
    # GROUP ( gene/funfam/gene-funfam all under col SPECIES with col SPECIES_TYPE saying what refers to...)                            
    if ( analysis_type == 'timing' ) { 
        # TIMING ONLY
        df_analysis <- group_by( df_analysis, 
                                    SPECIES_TYPE, 
                                    DATA_SOURCE, 
                                    TXP, 
                                    CANCER_TYPE, 
                                    OVERALL_TIMING, 
                                    SPECIES ) 
        # Remove cols not req in select
        df_analysis <- select( df_analysis, -c( ANY_REGION_DUPLICATED, PYCLONECLONAL_SC ) )
    
    } else if ( analysis_type == 'timing-clonality' ) {
        # TIMING & CLONALITY (really just used for "clonal only" plots)
        # remove any NA for clonality as these will not work well with diversity plots
        df_analysis <- df_analysis %>%
            filter( !is.na( PYCLONECLONAL_SC ) )
        df_analysis <- group_by( df_analysis, 
                                    SPECIES_TYPE, 
                                    DATA_SOURCE, 
                                    TXP, 
                                    CANCER_TYPE, 
                                    OVERALL_TIMING, 
                                    PYCLONECLONAL_SC, 
                                    SPECIES )
        # Remove cols not req in select
        df_analysis <- select( df_analysis, -c( ANY_REGION_DUPLICATED ) )

    } else if ( analysis_type == 'duplicated-clonality' ) {
        # ANY REGION DUPLICATED & CLONALITY
        # remove any NA as these will not work well with diversity plots
        df_analysis <- df_analysis %>%
            filter( !is.na( ANY_REGION_DUPLICATED ) )
        df_analysis <- group_by( df_analysis, 
                                    SPECIES_TYPE, 
                                    DATA_SOURCE, 
                                    TXP, 
                                    CANCER_TYPE, 
                                    ANY_REGION_DUPLICATED, 
                                    PYCLONECLONAL_SC, 
                                    SPECIES )
        # Remove cols not req in select
        df_analysis <- select( df_analysis, -c( OVERALL_TIMING ) )
    } else {
        print( "analysis_type should be 'timing', 'timing-clonality' or 'duplicated-clonality'" )
        return( -1 )
    }

    # total muts and distinct
    df_analysis <- df_analysis %>%
        mutate( tot_muts = sum( NUM_MUTS , na.rm = TRUE ) ) %>% 
        select( -c( NUM_MUTS ) ) %>%
        distinct( )

    return( df_analysis )
}

# Pivot wider so have comparisons (early/late etc) in cols for each species
# analysis_type: 'timing', 'timing-clonality', 'duplicated-clonality'
pivot_analysis <- function( dfa, analysis_type = 'timing', names_separator = '_' ) {
    if ( analysis_type == 'timing' ){
        df_pivot <- pivot_wider( dfa, 
                                    names_from = c( CANCER_TYPE, SPECIES_TYPE, OVERALL_TIMING ), 
                                    names_sep = names_separator,
                                    values_from = tot_muts, 
                                    values_fill = 0 )
    } else if ( analysis_type == 'duplicated-clonality' ) {
        df_pivot <- pivot_wider( dfa, 
                                    names_from = c(  CANCER_TYPE, SPECIES_TYPE, ANY_REGION_DUPLICATED, PYCLONECLONAL_SC ), 
                                    names_sep = names_separator,
                                    values_from = tot_muts, 
                                    values_fill = 0 )
    } else if ( analysis_type == 'timing-clonality' ) {
        df_pivot <- pivot_wider( dfa, 
                                    names_from = c(  CANCER_TYPE, SPECIES_TYPE, OVERALL_TIMING, PYCLONECLONAL_SC ), 
                                    names_sep = names_separator,
                                    values_from = tot_muts, 
                                    values_fill = 0 )
    } else {
        print( "analysis_type should be 'timing', 'timing-clonality', 'duplicated-clonality' " )
        return( -1 )
    }
    return(df_pivot)
}

# convert to suitable iNEXT input format
convert_for_iNEXT <- function( dfa ){
    dfi <- ungroup( dfa )
    dfi <- select( dfi, -c( DATA_SOURCE, TXP ) )
    dfi <- column_to_rownames( dfi, var = 'SPECIES' ) 
    
    # return data.frame (column_to_rownames coverts tibbles to df implicitly)
    return( dfi )
}

# Convert iNext df objects to log form for the diversity and confidence intervals
# assemblage = 'size' or 'coverage' for type of diversity table to convert
# dual_axis_log = FALSE : also convert the x axis (num species / m)
convert_inext_to_log <- function( df_inext, dual_axis_log = FALSE, log_base = 'natural' ){
    # convert both size and coverage assemblage 
    
    # size
    ln_inext <- conv_log( df_inext$iNextEst$size_based, dual_axis_log, log_base )
    if ( ln_inext != -1 ){
        df_inext$iNextEst$size_based <- ln_inext
    }
    
    # coverage
    ln_inext <- conv_log( df_inext$iNextEst$coverage_based, dual_axis_log, log_base )
    if ( ln_inext != -1 ) {
        df_inext$iNextEst$coverage_based <- ln_inext
    }

    # AsyEst
    ln_asy <- conv_asy( df_inext$AsyEst, log_base )
    if ( ln_asy != -1 ) {
        df_inext$AsyEst <- ln_asy
    }
    
    return( df_inext )
}

conv_log <- function( df_inext, dual_axis_log = FALSE, log_base = 'natural' ){
    # "private" function called from convert_inext_to_log
    
    # log convert
    ln_inext <- df_inext
    if ( log_base == 'log10' ){
        # log10
        ln_inext <- ln_inext %>%
            mutate( qD = log10(qD) ) %>%
            mutate( qD.LCL = log10(qD.LCL) ) %>%
            mutate( qD.UCL = log10(qD.UCL) )
        # also convert 'x' axis?
        if ( dual_axis_log == TRUE ){
            ln_inext <- mutate( ln_inext, m = log10(m) )
        }
    }else if ( log_base == 'log2' ){
        # log2
        ln_inext <- ln_inext %>%
            mutate( qD = log2(qD) ) %>%
            mutate( qD.LCL = log2(qD.LCL) ) %>%
            mutate( qD.UCL = log2(qD.UCL) )
        # also convert 'x' axis?
        if ( dual_axis_log == TRUE ){
            ln_inext <- mutate( ln_inext, m = log2(m) )
        }
    }else if ( log_base == 'natural' ){
        # natural log, e
        ln_inext <- ln_inext %>%
            mutate( qD = log(qD) ) %>%
            mutate( qD.LCL = log(qD.LCL) ) %>%
            mutate( qD.UCL = log(qD.UCL) )
        # also convert 'x' axis?
        if ( dual_axis_log == TRUE ){
            ln_inext <- mutate( ln_inext, m = log(m) )
        }
    }else {
        print( "invalid log_base, use 'natural', 'log10' or 'log2' ")
        return( -1 )
    }
    return(ln_inext)
}

# just converts the AsyEst which has different var names 
# could be more elegant but didn't want to rename vars in case of confusion
conv_asy <- function( ln_asy, log_base = 'natural' ){
    
    if ( log_base == 'log10' ){
        # log10
        ln_asy <- ln_asy %>%
            mutate( Observed = log10(Observed) ) %>%
            mutate( Estimator = log10(Estimator) ) %>%
            mutate( LCL = log10(LCL) ) %>%
            mutate( UCL = log10(UCL) )
    }else if ( log_base == 'log2' ){
        # log2
        ln_asy <- ln_asy %>%
            mutate( Observed = log2(Observed) ) %>%
            mutate( Estimator = log2(Estimator) ) %>%
            mutate( LCL = log2(LCL) ) %>%
            mutate( UCL = log2(UCL) )
    }else if ( log_base == 'natural' ){
        # natural log, e
        ln_asy <- ln_asy %>%
            mutate( Observed = log(Observed) ) %>%
            mutate( Estimator = log(Estimator) ) %>%
            mutate( LCL = log(LCL) ) %>%
            mutate( UCL = log(UCL) )
   
    }else {
        print( "invalid log_base, use 'natural', 'log10' or 'log2' ")
        return( -1 )
    }
    return(ln_asy)
}

# ---------------------
# iNEXT plots
# ---------------------
inext_plot <- function( idst, atyp, data_file_type, an_dir, label_full, logplot=TRUE, log_type='log2', facetvar = "None" ) {
  # plot_types: three types of plots: sample-size-based rarefaction/extrapolation curve (type = 1); sample completeness curve (type = 2); coverage-based rarefaction/extrapolation curve (type = 3).
  plot_types <- c(1,2) # just use this value
  # atyp: 'timing', 'duplicated-clonality', 'timing-clonality'

  # Log plots?
  if ( logplot==TRUE ) {
    # If log2-log2 plots
    plot_info <- 'log2-log2-plot_'
    # log_type <- 'log2'
    plot_dst <- convert_inext_to_log( idst, dual_axis_log = TRUE, log_base = log_type )
  } else {
    # If not log plots (eg just FIE)
    plot_info <- 'plot_'
    log_type <- ''
    plot_dst <- idst
  }
  
  # Make plot!
  inext_plot = ggiNEXT3D( plot_dst, type=plot_types, facet.var=facetvar )

  # Customise axes
  if (logplot == TRUE){
    x_text <- paste0( log_type, '(number of mutations)' )
    y_text <- paste0( log_type, '(Hill-Shannon diversity) grouped by ', data_file_type )

  } else{
    x_text <- 'number of mutations' 
    y_text <- paste0( 'Hill-Shannon diversity grouped by ', data_file_type )

  }
  samp_re <- inext_plot[[1]]  # explicit samp size re
  samp_completeness <- inext_plot[[2]]  # explicit samp coverage
  
  # Customise axes
  # coverage/completeness
  samp_completeness <- samp_completeness +
              xlab( x_text ) +
              ylab( 'coverage proportion ' ) +
              theme( legend.position = "right" ) +
              ggtitle( paste0( paste0( cancers, collapse = '-' ), 
                                  " sample coverage ",
                                  " - ", str_replace_all( label_info, '-', ' ' ) )
                      )
  # diversity
  samp_re <- samp_re +
              xlab( x_text ) +
              ylab( y_text ) +
              theme( legend.position = "right" ) +
              ggtitle( paste0( paste0( cancers, collapse = '-' ), 
                                  " diversity by ", atyp,
                                  " - ", str_replace_all( label_info, '-', ' ' ) ),
                      subtitle = paste0( "Sample-sized based sampling curves (", paste0( conf * 100 ), "% CI)" ) 
                      ) 

  # if duplicated-clonal then customise shape as have 8 lines ( {backgr, FIE} x {duplicated, not} x {C,S} )
  # empty shapes background, filled shapes for FIEs; circles (non-dupe), triangles (duplicated)
  # handy guide: https://stackoverflow.com/questions/16813278/cycling-through-point-shapes-when-more-than-6-factor-levels  [firefox/notion: #ggplot #R]
  if ( atyp == "duplicated-clonality" ){
      samp_re <- samp_re +            
              scale_shape_manual( values = c( 20, 8, 20, 8, 20, 8, 20, 8  )  )   # simpler filled circle for clonal, mult-cross sub-clonal
              
      # just duplicated / not
      # http://www.cookbook-r.com/Graphs/Colors_(ggplot2) 
      cbPalette <- c( "#56B4E9", "#56B4E9", "#E69F00", "#E69F00",  "#0072B2","#0072B2", "#D55E00", "#D55E00" )  # blue/orange not/dupe (darker for FIEs)

      samp_re <- samp_re +            
          scale_colour_manual(values=cbPalette) +
          scale_fill_manual(values=cbPalette)
  }

  # save plots (1: R/E sample size, 2: samp completeness)
  ggsave( paste0( plot_info, label_full, '_sampsize_r-e.pdf' ), samp_re, path = an_dir )
  ggsave( paste0( plot_info, label_full, '_plot_sampcompleteness.pdf' ), samp_completeness, path = an_dir )
  
}  