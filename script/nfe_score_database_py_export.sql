-- data/nfe_score_database_py_export.sql
-- Nov 2023
-- for submission - this is used to export from db:
-- 		data/nfe_score_database_py_funvar_archive_nfe_pfh_mutfam_090_exported.tsv
-- post process:
-- cat nfe_score_database_py_funvar_archive_nfe_pfh_mutfam_090_exported.tsv |cut -f1-47,49-66 > nfe_score_database_py_funvar_archive_nfe_pfh_mutfam_090_exported_no_mut_id.tsv
-- rm data/nfe_score_database_py_funvar_archive_nfe_pfh_mutfam_090_exported.tsv

-- A. For CATH PDB FunFams
    SELECT DISTINCT
        nfe_type,
        nfe_version,
        data_source,
        TxP,
        cancer_type,
        variant_type,
        variant_class,
        vm_synonymous,
        sf_id,
        ff_id,
        rep_id,
        mutfam,
        mutfam_version,
        mfc_mut_count_missense,
        mfc_mut_count_silent,
        source_hugo_symbol,
        vm_uniprot_accession,
        vm_seq_no,
        vm_aa_change,
        pdb_res,
        on_scons_90,
        on_mcsa,
        on_ppi,
        on_nuc,
        on_lig,
        near_angstroms,
        near_scons_90,
        near_mcsa,
        near_nuc,
        near_lig,
        num_gd,
        num_ngd,
        num_na_gd,
        num_timing_early,
        num_timing_late,
        num_timing_subclonal,
        num_timing_unknown,
        num_timing_na,
        num_clonal,
        num_subclonal,
        num_naclonal,
        gnomad_af_full,
        gnomad_af_converted,
        diseases,
        cadd_mark,
        protein_cat_res,
        protein_disulphide,
        mutation_id,
        hid,
        mc_num_muts,
        num_patients,
        mut_count_sum_p,
        mut_count_sum_p_corr,
        mut_count_sum_p_corr_sig,
        weighted_mut_sum_p,
        weighted_mut_sum_p_corr,
        weighted_mut_sum_p_corr_sig,
        groupid,
        runid,
        taskid,
        mc_num_swissprot_ff,
        mc_num_cgc_genes,
        mc_num_res_ff_rep,
        mc_num_mut_res_ff_rep,
        mc_sig_mut_count_sum_corr,
        mc_sig_weighted_mut_sum_corr   
    FROM
        funvar_archive.nfe_pfh_mutfam_090
    WHERE
        vm_synonymous='FALSE' 
        AND ( ( nfe_type = 'MC' AND groupid = 'mc07' ) OR ( nfe_type='PFH' AND num_patients > 2 ) )


-- B. For AlphaFold / TED 
	-- 1) Created gnomad_af_converted
        ALTER TABLE funvar_archive.nfe_pfh_mutfam_090_af 
        ADD ( gnomad_af_converted NUMBER DEFAULT -1 );
        
        -- -- If Tx then extract main AF from gnomAD_AF=6.696e-06,gnomAD_NFE_AF=1.763e-05 (TCGA VarMap was different)
        -- SELECT mutation_id, 
        --     SUBSTR( gnomad_af_full, 11, INSTR( gnomad_af_full,',gnomAD' ) - 11 ) gnomad_af_converted,
        --     gnomad_af_full, gnomad_af_converted
        -- FROM funvar_archive.nfe_pfh_mutfam_090_af
        -- WHERE data_source='Tx';        

        UPDATE funvar_archive.nfe_pfh_mutfam_090_af
        SET gnomad_af_converted = SUBSTR( gnomad_af_full, 11, INSTR( gnomad_af_full,',gnomAD' ) - 11 )
        WHERE data_source = 'Tx' AND gnomad_af_full <>'-';
		-- 44 rows updated.

        -- TCGA just copy unless '-'
        UPDATE funvar_archive.nfe_pfh_mutfam_090_af
        SET gnomad_af_converted = gnomad_af_full
        WHERE data_source='TCGA' AND gnomad_af_full <>'-';
        -- 713 rows updated.

    -- 2) export date for FIE scoring [script/fie_scoring/nfe_main.py]
    SELECT
        nfe_type,
        nfe_version,
        data_source,
        txp,
        cancer_type,
        variant_type,
        variant_class,
        vm_synonymous,
        sf_id,
        ff_id,
        rep_id,
        mutfam,
        mutfam_version,
        mfc_mut_count_missense,
        mfc_mut_count_silent,
        source_hugo_symbol,
        vm_uniprot_accession,
        vm_seq_no,
        vm_aa_change,
        pdb_res,
        on_scons_90,
        on_mcsa,
        on_ppi,
        on_nuc,
        on_lig,
        on_lig_id,
        near_angstroms,
        near_scons_90,
        near_mcsa,
        near_nuc,
        near_lig,
        near_lig_id,
        num_gd,
        num_ngd,
        num_na_gd,
        num_timing_early,
        num_timing_late,
        num_timing_subclonal,
        num_timing_unknown,
        num_timing_na,
        num_timing_early_mphase,
        num_timing_late_mphase,
        num_timing_subclonal_mphase,
        num_timing_unknown_mphase,
        num_timing_na_mphase,
        num_clonal,
        num_subclonal,
        num_naclonal,
        num_clonal_mphase,
        num_subclonal_mphase,
        num_naclonal_mphase,
        gnomad_af_full,
        gnomad_af_converted,
        diseases,
        cadd_mark,
        protein_cat_res,
        protein_disulphide,
        mutation_id,
        hid,
        mc_num_muts,
        mut_count_sum_p,
        mut_count_sum_p_corr,
        mut_count_sum_p_corr_sig,
        weighted_mut_sum_p,
        weighted_mut_sum_p_corr,
        weighted_mut_sum_p_corr_sig,
        groupid,
        runid,
        taskid,
        mc_num_swissprot_ff,
        mc_num_cgc_genes,
        mc_num_res_ff_rep,
        mc_num_mut_res_ff_rep,
        mc_sig_mut_count_sum_corr,
        mc_sig_weighted_mut_sum_corr,
        num_muts,
        num_tumours,
        num_patients
    FROM
	    funvar_archive.nfe_pfh_mutfam_090_af
    WHERE
        vm_synonymous='FALSE' 
        AND ( 
                ( nfe_type = 'MC' AND groupid = 'mcaf02_03_04' ) 
            OR  ( nfe_type = 'PFH' AND num_tumours > 2 ) 
        )
