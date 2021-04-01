# Introduction

## Background 
Genotyping array is the most economical approach for conducting large-scale genome-wide genetic association studies. Thorough quality control is key to generate high integrity genotyping data and robust results. Quality control of genotyping array is generally a complicated process, as it requires intensive manual labor in implementing the established protocols and curating a comprehensive quality report. There is an urgent need to reduce manual intervention via an automated quality control process. 

## Results

Based on previously established protocols and strategies, we developed an R package GTQC (GenoTyping Quality Control), which manages to automate a majority of the quality control steps and produce a detailed HTML report comprising tables and figures. Here, we describe the concepts underpinning GTQC and demonstrate its effectiveness using a real genotyping dataset.

## Conclusion
R package GTQC streamlines a majority of the quality control steps and produces a detailed HTML report on a plethora of quality control metrics, thus enabling a swift and rigorous data quality inspection prior to downstream GWAS and related analyses. By significantly cutting down on the time on genotyping quality control procedures, GTQC ensures maximum utilization of available resources and minimizes waste and inefficiency.     


# Example

## Example Data 1: 1000 Genome data

1000G plink data downloaded from https://www.cog-genomics.org/plink/1.9/resources . 
Race information downloaded from https://www.internationalgenome.org/faq/can-i-get-phenotype-gender-and-family-relationship-information-samples/
