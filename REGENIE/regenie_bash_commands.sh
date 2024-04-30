export MKL_DEBUG_CPU_TYPE=5
./regenie \
--step 1 \
--bed /Users/arjhunswaminathan/Documents/N27895_M516348_C2_P6_B126/REGENIE/output \
--covarFile /Users/arjhunswaminathan/Documents/N27895_M516348_C2_P6_B126/REGENIE/covariates.txt \
--phenoFile /Users/arjhunswaminathan/Documents/N27895_M516348_C2_P6_B126/REGENIE/phenotype.txt \
--bsize 4098 \
--force-qt --lowmem \
--cv 5 \
--lowmem-prefix tmp_rg \
--ignore-pred \
--out /Users/arjhunswaminathan/Documents/N27895_M516348_C2_P6_B126/REGENIE/fit_bin_out

./regenie \
--step 2 \
--bgen /Users/arjhunswaminathan/Documents/N27895_M516348_C2_P6_B126/REGENIE/output.bgen \
--covarFile /Users/arjhunswaminathan/Documents/N27895_M516348_C2_P6_B126/REGENIE/covariates.txt \
--phenoFile /Users/arjhunswaminathan/Documents/N27895_M516348_C2_P6_B126/REGENIE/phenotype.txt \

--bsize 4098 --force-qt \
--cv 5 \
--pred /Users/arjhunswaminathan/Documents/N27895_M516348_C2_P6_B126/REGENIE/fit_bin_out_pred.list \
--ignore-pred \
--out /Users/arjhunswaminathan/Documents/N27895_M516348_C2_P6_B126/REGENIE/test_bin_out_step_2
