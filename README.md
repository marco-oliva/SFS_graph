# Sample Specific Strings Graph

Computes the Sample Specific Strings with `moni` and feed them into `themisto`.

To run:

```bash 
git clone https://github.com/marco-oliva/SFS_graph
cd SFS_graph
sudo docker build . -t sfs_container
sudo docker run -v <your data folder>:/data --rm -it --entrypoint bash sfs_container

# Now from the container
SFS_graph.py -r 19.fa -p selected_HG002SFS.asm.ec.fa
```
