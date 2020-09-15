MToolBox-variant-calling
========================

This wrapper performs QC, quality trimming of raw reads, read alignment, alignment filtering, variant calling. The final output is a VCF file.

What should I look to, here?
----------------------------

**TL;DR**

- the VCF file in the :code:`results/vcf` folder
- the BED file(s) in the :code:`results/<sample>` folder.

Both file formats can be imported in a genome browser (*eg* IGV) to visually inspect your results.

The VCF output
^^^^^^^^^^^^^^

The VCF (variant call format) file is roughly a table where, after a ton of comment lines (starting with :code:`##`), rows are variants and columns are samples. You can find a detailed description of the VCF format `here`_, although MToolBox-snakemake provides a slightly different version to report the allele heteroplasmy frequency. To spare you some headache, we'll give you a brief summary of the genotype info you'll find for each sample.

The :code:`FORMAT` field lists data types and order that are available for samples (following fields). Each sample has colon-separated data corresponding to the types specified in the :code:`FORMAT`.

- :code:`GT` (genotype) reports all the alleles found for that sample, where :code:`0` is the reference allele (:code:`REF` field) and :code:`1`, :code:`2`, ... are the alleles in the same order as in the :code:`ALT` field.
- :code:`DP` (depth) reports the total coverage depth for that site in the genome, *i.e.* the number of reads mapping on it.
- :code:`HF` (heteroplasmic frequency) reports, for each allele in :code:`GT` (excluding :code:`REF`), the HF.
- :code:`CILOW` (confidence interval, lower bound) reports, for each allele in :code:`GT` (excluding :code:`REF`), the lower bound of the CI.
- :code:`CIUP` (confidence interval, upper bound) reports, for each allele in :code:`GT` (excluding :code:`REF`), the upper bound of the CI.
- :code:`SDP` (strand read depth) reports, for each allele in :code:`GT` (excluding :code:`REF`), the number of times the allele was observed on the plus strand and on the minus strand, semi-colon separated.

If you consider this example:

.. code-block:: bash
    
    #CHROM       POS    ID  REF  ALT  QUAL  FILTER  INFO       FORMAT                   Scer_mt_500K                     Scer_mt_100K
    NC_001224.1  50000  .   A    C    .     PASS    AN=4;AC=2  GT:DP:HF:CILOW:CIUP:SDP  0/1:359:0.46:0.409:0.511:90;75   0/1:75:0.387:0.284:0.5:11;18
    NC_001224.1  69045  .   C    T    .     PASS    AN=2;AC=1  GT:DP:HF:CILOW:CIUP:SDP  0/1:365:0.033:0.018:0.057:0;12   ./.:.:.:.:.:.

sample Scer_mt_500K, in mt position 50000, has 359 aligned reads and one variant allele (:code:`C`) with HF=0.46. The CI for this HF is 0.409 to 0.511. The variant allele is supported by 90 reads on the plus strand and 75 reads on the minus strand.

The BED output
^^^^^^^^^^^^^^

The `BED (browser extensible data) file`_ is a useful and intuitive way to inspect the variant calling results through a genome browser. Once you import this file in a genome browser, variants will be colour-coded (blue for mutations, green for insertions, red for deletions) and shaded according to the HF (the darker the shade, the higher the HF).

What to do next?
^^^^^^^^^^^^^^^^

Once you have run the wrapper, you will notice that the :code:`results` folder includes a ton of files. A guide to these files is coming soon.

.. _`here`: https://www.internationalgenome.org/wiki/Analysis/vcf4.0
.. _`BED (browser extensible data) file`: https://m.ensembl.org/info/website/upload/bed.html