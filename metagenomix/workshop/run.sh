WORKDIR="C:\folder\where\to\run\metagenomix"

metagenomix create \
	--fastq-dir-illumina $WORKDIR/fastqs \
	--output-dir $WORKDIR/outputs \
	--metadata $WORKDIR/metadata.tsv \
	--pipeline $WORKDIR/configs/pipeline.txt \
	--databases $WORKDIR/configs/databases.yml \
	--co-assembly $WORKDIR/configs/coassembly.yml \
	--user-params $WORKDIR/configs/user_params.yml \
	--modules $WORKDIR/configs/modules.yml \
	--project-name my_project_name \
	--chunks 4 \
	--account nn9745k \
	--scratch \
	--no-jobs \
	--dev

