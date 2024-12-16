export MKL_DEBUG_CPU_TYPE=5
./regenie \
--step 1 \
--bed ../test_site/N27895_M516348_C2_P6_B126/REGENIE/output \
--covarFile ../test_site/N27895_M516348_C2_P6_B126/REGENIE/covariates.txt \
--phenoFile ../test_site/N27895_M516348_C2_P6_B126/REGENIE/phenotype.txt \
--bsize 4098 \
--force-qt --lowmem \
--cv 5 \
--lowmem-prefix tmp_rg \
--ignore-pred \
--out ../test_site/N27895_M516348_C2_P6_B126/REGENIE/fit_bin_out

./regenie \
--step 2 \
--bgen ../test_site/N27895_M516348_C2_P6_B126/REGENIE/output.bgen \
--covarFile ../test_site/N27895_M516348_C2_P6_B126/REGENIE/covariates.txt \
--phenoFile ../test_site/N27895_M516348_C2_P6_B126/REGENIE/phenotype.txt \

--bsize 4098 --force-qt \
--cv 5 \
--pred ../test_site/N27895_M516348_C2_P6_B126/REGENIE/fit_bin_out_pred.list \
--ignore-pred \
--out ../test_site/N27895_M516348_C2_P6_B126/REGENIE/test_bin_out_step_2
