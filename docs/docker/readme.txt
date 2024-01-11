# build an image
docker build -t gene_mining .

# run a container
docker run -dit --name gene_mining_teste -v /media/lorrana/LorranaHD/Github/domain-based-gene-mining:/domain_based_gene_mining gene_mining
