SigLoops:Fast, accuracy, less-memory and user-friendly loop caller
=======


![image](https://user-images.githubusercontent.com/57889560/113148126-11ec0f00-9264-11eb-8fda-6d0020aaf1a6.png)

Loops in Hi-C contact matrix showed enrichment in reads count over the background. The background 
of contact matrix followed a monotonically decreasing pattern along the vertical diagonal which was a 
function of the distance to the site of interest.

Usage
-----
    SigLoops usage: HiC loop calling first step: map pets to fragment space.
         -i          The input pets file dir.Pets file name format:"chr1.vs.chr2.pets".
                     Pets format:chr1, pos1, chr2, pos2.
         -f          The REsite fragment bed file.
         -o          The output dir.
         -c          The chromosome list file path to select to process chrs.
         -w          The fragment to make window. Default: 50.
         -s          The step fragment of window.Default:5.
         -p          The cpu used. Default:1. Warning: no more than5!
