import os

configfile: "config.yaml"

WORKDIR = os.getcwd()
ASSEMBLIES = config["ASSEMBLIES"]
DIAMOND_DB = config["DIAMOND_DB"]
DIAMOND_DB_DIR = config["DIAMOND_DB_DIR_PARENT"]
ID_U = config["ID_U"]

rule all:
	input:
		# expand("results_{assembly}/prokka_initial/{assembly}_sprok.tbl",
				#  assembly=ASSEMBLIES)
		expand("results_{assembly}/analysis/pos_corrected_debug.bed",
					assembly=ASSEMBLIES),
		expand("results_{assembly}/analysis/pos_initial_debug.bed",
					assembly=ASSEMBLIES)

rule prokka_annotation:
	input:
		"data/{assembly}.fasta"
	params:
		outdir="results_{assembly}/prokka_initial",
		PREFIX="{assembly}_sprok",
		pwd=WORKDIR,
		id_user=ID_U
	output:
		"results_{assembly}/prokka_initial/{assembly}_sprok.tbl",
		"results_{assembly}/prokka_initial/{assembly}_sprok.fsa"
	shell:
		"docker run --rm --user {params.id_user}:{params.id_user} -v {params.pwd}:/db staphb/prokka:latest prokka \
		--force --kingdom Bacteria \
		--outdir /db/{params.outdir} --prefix {params.PREFIX} \
		--genus Pantoea --locustag LZP \
		/db/{input}"

rule index_genome:
	input:
		"results_{assembly}/prokka_initial/{assembly}_sprok.fsa"
	output:
		"results_{assembly}/prokka_initial/{assembly}_sprok.fsa.fai"
	shell:
		"samtools faidx {input}"


rule get_nucl_genes_from_tbl:
	input:
		genome_mock_idx="results_{assembly}/prokka_initial/{assembly}_sprok.fsa.fai",
		genome="results_{assembly}/prokka_initial/{assembly}_sprok.fsa",
		sample="results_{assembly}/prokka_initial/{assembly}_sprok.tbl"
	output:
		"results_{assembly}/prokka_initial/{assembly}_sprok_genes.fa"
	shell:
		"./target/debug/tbltk tofasta \
		--genome {input.genome} {input.sample} > {output}"


rule diamond_search:
	input:
		genes="results_{assembly}/prokka_initial/{assembly}_sprok_genes.fa"
	benchmark:
		"results_{assembly}/prokka_initial/benchmark_diamond.txt"
	output:
		"results_{assembly}/prokka_initial/output_diamond.txt"
	params:
		diamonddb=DIAMOND_DB,
		id_user=ID_U,
		diamonddb_dir=DIAMOND_DB_DIR
	threads:
		12
	shell:
		"docker run --user {params.id_user}:{params.id_user} -v {params.diamonddb_dir}:/dbnr -v $(pwd):/db --rm diamond blastx -d /dbnr/{params.diamonddb} \
		-q /db/{input.genes} -F 15 -p {threads} \
		--matrix PAM30 -e 1e-20 -k 1 -o /db/{output} \
		--fast --outfmt 6 qseqid stitle sseqid qstart qend sstart send qframe btop"

rule fixing_frameshift:
	input:
		genome="results_{assembly}/prokka_initial/{assembly}_sprok.fsa",
		diamond="results_{assembly}/prokka_initial/output_diamond.txt",
		genes="results_{assembly}/prokka_initial/{assembly}_sprok_genes.fa",
		tbl="results_{assembly}/prokka_initial/{assembly}_sprok.tbl"
	benchmark:
		"results_{assembly}/framerust_res/benchmark.txt"
	threads:
		8
	output:
		cor="results_{assembly}/framerust_res/{assembly}_fixed.fa",
		debug_seq="results_{assembly}/framerust_res/debug.fa"
	shell:
		"./target/debug/uniref --debug -a {input.genome} \
		-b {input.diamond} -g {input.genes} -o {output.cor} \
		{input.tbl} -t {threads}"

rule prokka_new_annotation:
	input:
		"results_{assembly}/framerust_res/{assembly}_fixed.fa"
	params:
		outdir="results_{assembly}/prokka_corrected",
		PREFIX="{assembly}_eprok",
		pwd=WORKDIR,
		id_user=ID_U
	output:
		"results_{assembly}/prokka_corrected/{assembly}_eprok.tbl",
		"results_{assembly}/prokka_corrected/{assembly}_eprok.fsa"
	shell:
		"docker run --rm --user {params.id_user}:{params.id_user} -v {params.pwd}:/db staphb/prokka:latest prokka \
		--force --kingdom Bacteria \
		--outdir /db/{params.outdir} --prefix {params.PREFIX} \
		--genus Pantoea --locustag LZP \
		/db/{input}"

rule index_genome_corrected:
	input:
		"results_{assembly}/prokka_corrected/{assembly}_eprok.fsa"
	output:
		"results_{assembly}/prokka_corrected/{assembly}_eprok.fsa.fai"
	shell:
		"samtools faidx {input}"

rule mapping_debug_initial:
	input:
		genome="results_{assembly}/prokka_initial/{assembly}_sprok.fsa",
		genome_fai="results_{assembly}/prokka_initial/{assembly}_sprok.fsa.fai",
		debug_seq="results_{assembly}/framerust_res/debug.fa"
	output:
		"results_{assembly}/analysis/pos_initial_debug.bed"
	shell:
		"""
		minimap2 -cx asm5 --cs=long -t 1 {input.genome} \
			{input.debug_seq} | awk "{{ if( !(\$12<60) ) print \$0 }}" | paftools.js splice2bed - | \
			cut -f1,2,3,4,5,6  > {output}
		"""

rule mapping_debug_corrected:
	input:
		genome="results_{assembly}/prokka_corrected/{assembly}_eprok.fsa",
		genome_fai="results_{assembly}/prokka_corrected/{assembly}_eprok.fsa.fai",
		debug_seq="results_{assembly}/framerust_res/debug.fa"
	output:
		"results_{assembly}/analysis/pos_corrected_debug.bed"
	shell:
		"""
		minimap2 -cx asm5 --cs=long -t 1 {input.genome} \
		{input.debug_seq} | "{{ if( !(\$12<60) ) print \$0 }}" | paftools.js splice2bed - | \
		cut -f1,2,3,4,5,6  > {output}"
		"""

