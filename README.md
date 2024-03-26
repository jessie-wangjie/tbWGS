# tbWGS

A Nextflow pipeline for post-processing WGS data.
- DRAGEN Germline Calling
- Download the result from Basespace
  
The pipeline:
1. Annotate small variants from *.hard-filtered.vcf.gz using snpEff and output the variants with high risk for given target regions.
3. Annotate CNV variants from *.cnv.vcf.gz using snpEff and output the variants with high risk for given target regions.
4. Annotate  variants from *.sv.vcf.gz using snpEff and output the variants with high risk for given target regions.

## Running the pipeline
Run the pipeline with:

```console
nextflow run . --dragen_path <Path download from basespace> --target_bed <bed file>
```
