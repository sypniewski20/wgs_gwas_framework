Channel
	.fromPath("${params.input}")
	.set { ch_input }


process vcf_to_plink {
  publishDir "${baseDir}/${params.project_name}/temp"
	input:
	file(vcf) from ch_input
	output:
	tuple val(params.project_name),
	file("${params.project_name}_maf.bed"),
	file("${params.project_name}_maf.bim"),
	file("${params.project_name}_maf.fam"),
  file("${params.project_name}_maf.log") into plink_format
	script:
	"""
  plink --vcf ${vcf} --maf 0.05 \\
  --geno 0.05 \\
  --vcf-half-call r \\
  --make-bed \\
  --out ${params.project_name}_maf \\
  --const-fid 0 \\
  --threads ${params.threads} \\
  --hwe 10e-06 \\
  --${params.species} \\
  --update-sex ${baseDir}/${params.sex_info} \\
  --${params.include_chr}
	"""
}

plink_format
  .into { ld_est; filtering }



process ld_est {
  publishDir "${baseDir}/${params.project_name}/temp"
	input:
	file(plink) from ld_est
	output:
	tuple val(params.project_name),
	file("${params.project_name}_ld.prune.in"),
  file("${params.project_name}_ld.log") into prune_in
	script:
	"""
  plink --bfile ${params.project_name}_maf \\
  --indep 50 5 2 \\
  --${params.species} \\
  --out ${params.project_name}_ld
	"""
}

filtering
  .join(prune_in)
  .set{pruning}


process ld_pruning {
  publishDir "${baseDir}/${params.project_name}/temp"
	input:
	file(plink) from pruning
	output:
	tuple val(params.project_name),
	file("${params.project_name}_maf_pruned.bed"),
	file("${params.project_name}_maf_pruned.bim"),
	file("${params.project_name}_maf_pruned.fam"),
  file("${params.project_name}_maf_pruned.eigenvec"),
  file("${params.project_name}_maf_pruned.log") into pruned_plink
	script:
	"""
  plink --bfile ${params.project_name}_maf \\
  --threads ${params.threads} \\
  --extract ${params.project_name}_ld.prune.in \\
  --${params.species} \\
  --make-bed \\
  --out ${params.project_name}_maf_pruned \\
  --pca 20
	"""
}

process run_assoc {
  publishDir "${baseDir}/${params.project_name}/output"
	input:
	file(plink) from pruned_plink
	output:
	tuple val(params.project_name),
	file("*") into results
	script:
	"""
  plink --bfile ${params.project_name}_maf_pruned \\
  --logistic hide-covar sex perm \\
  --threads ${params.threads} \\
  --ci 0.95 \\
  --adjust qq-plot \\
  --pheno ${baseDir}/${params.pheno} \\
  --covar ${baseDir}/${params.covar} \\
  --covar-name ${params.covar_names} \\
  --${params.species} \\
  --1 \\
  --out ${params.project_name}_maf_pruned
	"""
}