## file manipulation of the trinity pipeline results

cut -f1-3,13 myTrinotate.tsv >> go.tsv#select the columns with only , transcript id, gene id, sp prot results and blast results
cut -f2,4 go.tsv >> go1.tsv #selecting onlt the tranxript id and the blast results
awk -F'\t' '$2 != "."' go1.tsv >go2.tsv## select only those transcripts with blstx results,i.e go terms


grep -Fwf ../diff go2.tsv > go3.txt ## find differentially expressed genes and their go terms
grep -Fwf ../diff go2.tsv > go4.txt #for deseq expressed genes: none were found
grep -Fwf biosynthsesisgenes.tsv go2.tsv > go5.txt #find the biosynthsesisgenes in the file with all goterms, all of them present
grep -Fwf biosynthsesisgenes.tsv go3.txt > go6.txt # find biosynthesis genes from Jens in the diff expressed genes, they are not there

grep -Fwf diffexpressedgenes1.tsv ../Thesis/Canuella_transcriptomics/results/go > go4.txt
