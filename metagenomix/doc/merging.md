:construction:

# Merging

This section will explain how to use the `metagenomix merge` command to 
combine the per sample or per co-assembly
[input units](https://github.com/FranckLejzerowicz/metagenomix/blob/main/metagenomix/doc/tutorials/code/io_unit.md) 
into:
* feature tables: "horizontal concatenation" of the input unit vectors, with 
  sample (or co-assembly group) names are columns. 
* long-format tables: "vertical concatenation" of the input unit outputs, 
  with the sample (or co-assembly group) added as a content to new column. 