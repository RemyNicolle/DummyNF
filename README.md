# DummyNF

To run the pipe, first modify the params in nextflow.config, in particular the following params:

- sampleInputDir: must be the path to the fq files
- sampleList: the path to the file listing the sample names and path
- outputDir: the path to a diretoy to be created with the output

Then run the pipe as folow:

```shell
~/nextflow  -c nextflow.config run .
```



The dummy pipe on c elegans takes about 2 CPU hours.





