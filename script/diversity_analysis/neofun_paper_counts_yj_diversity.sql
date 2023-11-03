-- git/funvar-tracerx/script/diversity_analysis/neofun_paper_counts_yj_diversity.sql
-- BASED ON script/v3.2/neofun_v3.2_paper_counts_yj_diversity.sql

-- split from _yj07 (especially part f)
-- neofun_v3.2_SI_yj
-- https://docs.google.com/spreadsheets/d/1nDk9dkMRZm5WJCZUP2PXn3f4P09e-khM1rkQMI11ssc/edit#gid=2032736972

-- Also: neofun/script/sql/vw/vw_diversity.sql
-- NOTE: GROUP BY queries with use SUM() to tally FIEs/muts as now based on view vw_diversity that does the COUNTs

-- A) Complete mut/FIE tables with all cols timing, clonal etc with species_type and species (ie gene/funfam/gene-funfam!)
    -- All mutations, inc variant_class, score_threshold and drivermut for filtering FIE/backgrou d etc
    -- No TxP filter, but only applies to Tx421 data anyway (i.e. overall_timing + other fields Tx-sspecific)
    --      data/diversity_yj07f_full_table_gene.csv
    --      data/diversity_yj07f_full_table_funfam.csv
    SELECT
        'gene' AS species_type,
        -- 'funfam' AS species_type,
        data_source,
        txp,
        variant_class,
        cancer_type,
        any_region_duplicated,
        pycloneclonal_sc,
        overall_timing,
        -- what is species type? gene/funfam
        source_hugo_symbol AS species,
        -- funfam AS species,
        score_threshold,
        drivermut,
        mutant_expressed,
        SUM(num_muts) num_muts,
        SUM(num_rows) num_rows
    FROM
        funvar_tx.vw_diversity
    -- filtering is done in R scripts
    -- WHERE 
    --     cancer_type IN ( 'LUAD', 'LUSC' )
    GROUP BY
        'gene',
        -- 'funfam',
        data_source,
        txp,
        variant_class,
        cancer_type,
        any_region_duplicated,
        pycloneclonal_sc,
        overall_timing,
        -- what is species type? gene/funfam
        source_hugo_symbol,
		-- funfam,
        score_threshold,
        drivermut,
        mutant_expressed
    ORDER BY
       'gene',
        -- 'funfam',
        data_source,
        txp,
        variant_class,
        cancer_type,
        any_region_duplicated,
        pycloneclonal_sc,
        overall_timing,
        -- what is species type? gene/funfam
        source_hugo_symbol,
		-- funfam,
        score_threshold,
        drivermut,
        mutant_expressed
    ;
