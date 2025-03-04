nextflow run workflows/WTS -preview -with-dag WTS.dot
dot -Tpng WTS.dot -o WTS-dag.png
