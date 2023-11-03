--------------------------------------------------------
--  File created - Wednesday-November-01-2023   
--------------------------------------------------------
--------------------------------------------------------
--  DDL for View VW_DIVERSITY
--------------------------------------------------------

  CREATE OR REPLACE FORCE EDITIONABLE VIEW "FUNVAR_TX"."VW_DIVERSITY" ("DATA_SOURCE", "TXP", "VARIANT_CLASS", "CANCER_TYPE", "ANY_REGION_DUPLICATED", "PYCLONECLONAL_SC", "OVERALL_TIMING", "SOURCE_HUGO_SYMBOL", "FUNFAM", "SCORE_THRESHOLD", "DRIVERMUT", "MUTANT_EXPRESSED", "NUM_MUTS", "NUM_ROWS") AS 
  SELECT
            fs.data_source,
            fs.txp,
            fs.variant_class,
            fs.cancer_type,
            ovt.any_region_duplicated,
            ovt.pycloneclonal_sc,
            ovt.overall_timing,
            fs.source_hugo_symbol,
            CONCAT( CONCAT( CONCAT( fs.sf_id, '_' ), fs.ff_id ), CONCAT( '_', fs.rep_id ) ) funfam,
            sig.score_threshold,
            dr.drivermut,
            ovt.mutant_expressed, -- from Dropbox/bioinf/neofun/data/in/alex/20230208 needs including on ovt!				
            COUNT( DISTINCT fs.mutation_id ) num_muts,  -- which can be FIEs depending on query type
            COUNT(*) num_rows
        FROM
            funvar_tx.mvw_tx_fs_3d fs
        -- Overall timing table Dec 22 (ITX_OVT03)
        INNER JOIN
            funvar_import.mut_tx_ov_timing ovt 
            ON fs.mutation_id = ovt.mutation_id
        -- Tx421 expressed genes (any expression level)		
        INNER JOIN (
            SELECT
                DISTINCT expressed_gene
            FROM
                funvar_import.expressed_genes_tx
            ) exp_tx
            ON exp_tx.expressed_gene = fs.source_hugo_symbol
        -- FIEs or Tx drivers
        -- use correct WHERE clause!
        -- MC3 FIEs (LEFT JOIN)
        LEFT JOIN
            (
                SELECT 
                    mutation_id, 
                    score_threshold 
                FROM 
                    funvar_tx.nfe_sig_mc_pfh 
                WHERE  score_threshold = 'MC3'
            ) sig
            ON sig.mutation_id = fs.mutation_id
        -- Tx driverMuts (LEFT JOIN)
        LEFT JOIN
            (
            SELECT 
                mutation_id,
                drivermut
            FROM
                funvar_import.mut_tx_drivers
            WHERE drivermut = 'TRUE'
            ) dr
            ON dr.mutation_id = fs.mutation_id
        WHERE 
            -- NOTE ONLY fs.pass_useforplot = 'TRUE' are correct FIEs - export includes all to see what missing
            fs.pass_useforplot = 'TRUE' -- (only correct ones: PASS in Tx or useforplot in TCGA-TxP)
            AND fs.data_source = 'Tx' 
            AND fs.overall_timing IN ( 'early', 'late', 'unknown' )
            -- 26 10 2022 important extra condition for Tx from inherited from yj05 
            AND NOT ( fs.txp='N' AND fs.overall_timing IS NULL )
        GROUP BY
            -- standard group for VIEW
            fs.data_source,
            fs.txp,
			fs.variant_class,
            fs.cancer_type,
            ovt.any_region_duplicated,
            ovt.pycloneclonal_sc,
            ovt.overall_timing,
            fs.source_hugo_symbol,
            CONCAT( CONCAT( CONCAT( fs.sf_id, '_' ), fs.ff_id ), CONCAT( '_', fs.rep_id ) ),
            sig.score_threshold,
            dr.drivermut,
            ovt.mutant_expressed
        ORDER BY
            -- standard ORDER for VIEW
            fs.data_source,
            fs.txp,
			fs.variant_class,
            fs.cancer_type,
            ovt.any_region_duplicated,
            ovt.pycloneclonal_sc,
            ovt.overall_timing,
            fs.source_hugo_symbol,
            CONCAT( CONCAT( CONCAT( fs.sf_id, '_' ), fs.ff_id ), CONCAT( '_', fs.rep_id ) ),
            sig.score_threshold,
            dr.drivermut,
            ovt.mutant_expressed
;
