# pyscenic

conda create -y -n pyscenic python=3.10
conda activate pyscenic

pip install numpy==1.23.4 # numpy verison problem
pip install pandas==1.2.4
pip install pyscenic


# download public motif_databases
destin="path_to_save_motif_databses"
cd $destin
wget https://resources.aertslab.org/cistarget/databases/homo_sapiens/hg38/refseq_r80/mc_v10_clust/gene_based/hg38_10kbp_up_10kbp_down_full_tx_v10_clust.genes_vs_motifs.rankings.feather ./
wget https://resources.aertslab.org/cistarget/databases/homo_sapiens/hg38/refseq_r80/mc_v10_clust/gene_based/hg38_10kbp_up_10kbp_down_full_tx_v10_clust.genes_vs_motifs.scores.feather ./
wget https://resources.aertslab.org/cistarget/databases/homo_sapiens/hg38/refseq_r80/mc_v10_clust/gene_based/hg38_500bp_up_100bp_down_full_tx_v10_clust.genes_vs_motifs.rankings.feather ./
wget https://resources.aertslab.org/cistarget/databases/homo_sapiens/hg38/refseq_r80/mc_v10_clust/gene_based/hg38_500bp_up_100bp_down_full_tx_v10_clust.genes_vs_motifs.scores.feather ./
wget https://resources.aertslab.org/cistarget/motif2tf/motifs-v10nr_clust-nr.hgnc-m0.001-o0.0.tbl ./
wget https://resources.aertslab.org/cistarget/tf_lists/allTFs_hg38.txt ./

# step1: GRN inference, generation of coexpression modules
conda activate pyscenic

public_data_dir="path_to_motif_databases"
input_dir="path_to_input"
output_dir="path_to_output"

for celltype_cur in cholangiocyte stellate_cell hepatocyte endothelial_cell;do
  echo $celltype_cur
  pyscenic grn ${input_dir}${celltype_cur}".loom" ${public_data_dir}"allTFs_hg38.txt" \
  -o ${output_dir}${celltype_cur}.adj.csv --num_workers 20

  pyscenic ctx ${output_dir}${celltype_cur}.adj.csv \
      ${public_data_dir}hg38_10kbp_up_10kbp_down_full_tx_v10_clust.genes_vs_motifs.rankings.feather \
      ${public_data_dir}hg38_500bp_up_100bp_down_full_tx_v10_clust.genes_vs_motifs.rankings.feather \
      --annotations_fname ${public_data_dir}motifs-v10nr_clust-nr.hgnc-m0.001-o0.0.tbl \
      --expression_mtx_fname ${input_dir}${celltype_cur}.loom \
      --output ${output_dir}${celltype_cur}.reg.csv \
      --mask_dropouts \
      --num_workers 20

  pyscenic aucell \
      ${input_dir}${celltype_cur}.loom \
      ${output_dir}${celltype_cur}.reg.csv \
      --output ${output_dir}${celltype_cur}_scenic_output.loom \
      --num_workers 20
done
