# LDSC
ldsc_dir="path_to_output"
plink_dir="path_to_chromosome_split_plink"
gwas_list=(Nat22_cALT_NAFLD.sumstats BBJ.Alb.sumstats BBJ.ALP.sumstats BBJ.ALT.sumstats BBJ.AST.sumstats BBJ.BUN.sumstats BBJ.Ca.sumstats BBJ.CK.sumstats BBJ.Cl.sumstats BBJ.CRP.sumstats BBJ.DBP.sumstats BBJ.GGT.sumstats BBJ.HbA1c.sumstats BBJ.HDL-C.sumstats BBJ.K.sumstats BBJ.LDL-C.sumstats BBJ.MAP.sumstats BBJ.Na.sumstats BBJ.Plt.sumstats BBJ.PP.sumstats BBJ.RBC.sumstats BBJ.sCr.sumstats BBJ.TBil.sumstats BBJ.TC.sumstats BBJ.TG.sumstats BBJ.TP.sumstats BBJ.UA.sumstats)
annot_list=('all_eqtl' "liver_eqtl" "ieqtl" "hepatocyte_all" 'hepatocyte_i' 'hepatocyte_liver' 'cholangiocyte_all' 'cholangiocyte_liver' 'cholangiocyte_i' 'stellate_cell_all' 'stellate_cell_liver' 'stellate_cell_i' 'endothelial_cell_all' 'endothelial_cell_liver' 'endothelial_cell_i' )

# make weight, frequency files that match or data
for chr in {1..22};do
  python ./ldsc.py \
  --bfile $plink_dir"scNAFLD_split."$chr \
  --l2 --ld-wind-cm 1 \
  --out $ldsc_dir"output/scNAFLD.weights."${chr}
done

for chr in {1..22};do
plink --bfile ${plink_dir}"scNAFLD_split."${chr} \
  --freq --keep-allele-order \
  --out ${ldsc_dir}"/output/scNAFLD_split."${chr}
done


# make annotation file, ld score
for annot in "${annot_list[@]}";do
  echo $annot
  if [ ! -d $ldsc_dir"output_"${annot} ];then
    mkdir $ldsc_dir"output_"${annot}
  fi

  for chr in {1..22};do
    echo $chr
    python make_annot.py \
    --bed-file $ldsc_dir"input/"${annot}"."${chr}".bed" \
    --bimfile $plink_dir"scNAFLD_split."$chr".bim" \
    --annot-file $ldsc_dir"output_"${annot}"/"${annot}"."${chr}".tmp.gz"

    zcat $ldsc_dir"output_"${annot}"/"${annot}"."${chr}".tmp.gz" | awk '{OFS="\t"; if (NR==1) {print "base",$1} else {print "1",$1} }' | bgzip -c > $ldsc_dir"output_"${annot}"/"${annot}"."${chr}".annot.gz"
    rm $ldsc_dir"output_"${annot}"/"${annot}"."${chr}".tmp.gz"
  done

  for chr in {1..22};do
    python ldsc.py \
    --l2 \
    --bfile $plink_dir"scNAFLD_split."$chr \
    --ld-wind-cm 1 \
    --annot ${ldsc_dir}"output_"${annot}"/"${annot}"."${chr}".annot.gz" \
    --thin-annot \
    --out $ldsc_dir"output_"${annot}"/"${annot}"."${chr}
  done
done


# heritability partitioning on various phenotype GWAS
gwas_dir='path_to_gwas_sumstats'

for gwas_cur in "${gwas_list[@]}";do
  echo $gwas_cur
  gwas_sumstats=${gwas_dir}${gwas_cur}

  for annot in "${annot_list[@]}";do
    echo $annot
    python ldsc.py \
      --h2 $gwas_sumstats \
      --ref-ld-chr $ldsc_dir"output_"${annot}"/"${annot}".@" \
      --out $ldsc_dir"output/"${gwas_cur}"__"${annot} \
      --w-ld-chr $ldsc_dir"output/scNAFLD.weights.@" \
      --frqfile-chr $ldsc_dir"output/scNAFLD_split.@" \
      --overlap-annot --print-cov --print-coefficients --print-delete-vals --not-M-5-50
  done
done
