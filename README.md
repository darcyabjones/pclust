# pclust

A pipeline to cluster and compare clusters of distantly related proteins.

pclust does lots of the heavy lifting involved in beginning to define protein families beginning with databases of millions of proteins.

Essentially it does the following:

1) Cluster all proteins into groups with 80% reciprocal coverage of each member.
   This means that most proteins should share the same domain structure,
   but may be highly sequence diverse.
2) Construct MSAs from the clustered protein sequences.
3) [Optionally] Construct rough trees for the clusters.
4) Construct an HHSuite database for each cluster using enriched MSAs.
5) Perform all-vs-all HMM-HMM comparisons of the clusters, to identify
   full-length remote homologues, truncated versions of proteins, and
   clusters that share common domains.
6) [Optionally] Search for HMM-HMM database matches to Pfam, SCOP, PDB, and Uniref/Uniclust
   for additional characterisation of clusters or matches between clusters.


Our focus in on defining relationships between groups of [Fungal Effector proteins](https://en.wikipedia.org/wiki/Effector_(biology))
which might have structural similarity but very low-sequence similarity.
Additional steps for analysing these relationships will be added in the future.

## Install

The pipeline is implemented in [Nextflow](https://www.nextflow.io/) and provides pre-built docker and singularity containers containing dependencies.

The pipeline uses:

- [Nextflow](https://www.nextflow.io/)
- [MMSeqs](https://github.com/soedinglab/MMseqs2)
- [HHSuite](https://github.com/soedinglab/hh-suite)
- [MAFFT](https://mafft.cbrc.jp/alignment/software/)
- [FastTree](http://www.microbesonline.org/fasttree/)
- [ffdb](https://github.com/darcyabjones/ffdb/)
- [ffindex](https://github.com/ahcm/ffindex) (Well actually we use the [HHSuite fork](https://github.com/soedinglab/hh-suite/tree/master/lib/ffindex), but credit to the original authors).


To run using the containers (recommended) you'll need to install [Singularity](https://sylabs.io/guides/latest/user-guide/) (recommended, no sudo required) or [Docker](https://docs.docker.com/install/) (sudo required).


Nextflow can automatically pull the appropriate container for you from docker/singularity-hub
if you provide the config option.
To pull the images manually:

```
singularity key pull 43791281018404BF71E1BAD5D520C428FF69B8C1
singularity pull library://darcyabjones/default/pclust:v0.0.1

# or

sudo docker pull darcyabjones/pclust:v0.0.1
```

Replacing `v0.0.1` with the version that you actually want.
Container versions are pinned to the pipeline releases and tags, so try to match them up.


If you want to squeeze every bit of performance out of your CPU, you can try to compile the software or containers yourself on your own hardware.
See the [containers](https://github.com/darcyabjones/pclust/tree/master/containers) folder for instructions on building containers.
We have made efforts to support both SSE and AVX2 instructions in the containers, so unless you've got some really fancy or really old gear the containers should do well.


## Running

Pclust is called using the `nextflow` command.
You can get a list of parameters like this.

```
# To have nextflow fetch the latest version from github.
nextflow run darcyabjones/pclust --help

# Download it yourself and execute from within the directory.
curl -sSL https://github.com/darcyabjones/pclust/archive/v0.0.1.tar.gz | tar -xf - -C .
cd pclust-v0.0.1
nextflow run ./main.nf --help
```

To run the pipeline using the containers (nextflow will pull the containers for you).

```
nextflow run -resume -profile singularity darcyabjones/pclust --help
nextflow run -resume -profile docker darcyabjones/pclust --help
```


Note for docker you may need to use `sudo` or modify the docker configuration to use `sudo`.
See [nextflow documentation](https://www.nextflow.io/docs/latest/config.html#config-docker) on docker configuration.



To run the full pipeline (without searching Pfam etc):

```
nextflow run \
  -resume \
  -profile singularity,nimbus \
  darcyabjones/pclust \
  --seqs my_proteins.fasta \
  --enrich_seqs uniref50.fasta \
  --tree
```


The `enrich_seqs` can be the same as the `seqs`, but this is usually only meaningful for when `seqs` already has a lot (millions) of diverse members (or for testing).
For clustering smaller datasets, try using an appropriate [Uniref](https://www.uniprot.org/help/uniref) enrichment database.


To run the searches against Pfam, Scop, PDB, and uniref you'll need to download those databases and (unfortunately) organise them into a particular folder structure.

You can download preformatted databases from <http://wwwuser.gwdg.de/~compbiol/data/hhsuite/databases/hhsuite_dbs/> and <https://uniclust.mmseqs.com/>.
Note for SCOP, make sure you download the one that says "hhsuite3" in the name (others are just clustered versions of the database not hmms).

Once downloaded and extracted, you can move or symlink the files to a structure that pclust expects.
The reason that we require a specific structure is because the databases are composed of 6 files, so it's easier to provide a folder containing the database than six individual options, and the databases have version information in the filenames, which makes hard-coding names difficult.

The folder should have the following structure:

```
hhdatabases/pfam_a3m.ffdata
hhdatabases/pfam_a3m.ffindex
hhdatabases/pfam_hhr.ffdata
hhdatabases/pfam_hhr.ffindex
hhdatabases/pfam_cs219.ffdata
hhdatabases/pfam_cs219.ffindex

hhdatabases/scop_a3m.ffdata
hhdatabases/scop_a3m.ffindex
            ...

hhdatabases/pdb_a3m.ffdata
hhdatabases/pdb_a3m.ffindex
            ...

hhdatabases/uniref_a3m.ffdata
hhdatabases/uniref_a3m.ffindex
            ...
```

The different databases (`pfam_`, `scop_`, `pdb_`, `uniref_`) don't have to be in the same folder (but it's somewhat convenient) and some can be excluded.
The name of the containing folder (here `hhdatabases`) is also not important, it can be anything.

Something like the following bash loop(s) will rename the files appropriately.

```
for f in scop70_1.75_*.ff{data,index}; do mv $f scop_${f##*_}; done
for f in pdb70_*.ff{data,index}; do mv $f pdb_${f##*_}; done
for f in uniclust30_2018_08_*.ff{data,index}; do mv $f uniref_${f##*_}; done
```


To run the HMM searches against these databases call pclust with the `--hhpfam`, `--hhscop`, `--hhpdb`, and `--hhuniref` options:

```
nextflow run \
  -resume \
  -profile singularity,nimbus \
  darcyabjones/pclust \
  --seqs my_proteins.fasta \
  --enrich_db uniref50.fasta \
  --hhpfam hhdatabases \
  --hhscop hhdatabases \
  --hhpdb hhdatabases
```

Would run searches against pfam, pdb and scop, but not uniref/uniclust even if it was present in the folder.


You can stop the pipeline at various points in the pipeline and provide precomputed options to skip sections.
Generally, it will be easier if you run all clustering steps within the pipeline.
But it is probably a good idea to run the hmm database searches on a compute cluster using FFindex
if your cluster doesn't have good support for running things like nextflow.


## Arguments

```
Mandatory Arguments:
  --seqs              The protein sequences to cluster.
  --db                The proteins to cluster as an MMseqs formatted database.
                      Either `--db` or `--seqs` must be provided to run clustering.
  --enrich_seqs       Sequences to enrich the clusters with for clustering
                      and converting into HHsuite databases.
  --enrich_db         A database to enrich clusters and hhsuite database with.
                      Should be an MMseqs formatted sequence database.
                      If neither `--enrich_seqs` or `--enrich_db` then the seqs
                      will be used for enrichment.

Options:
  --nomsa             Stop the pipeline after clustering.
  --tree              Predict quick phylogenetic trees for each MSA.
  --noremote          Stop the pipeline after MSAs, and don't run the HMM-HMM comparisons.
  --clusters          Provide existing clusters as an MMseqs database instead of using
                      the built-in pipeline. Note that the seq database must also be provided.
  --msas              Provide existing MSAs as an MMSeqs MSA database (fasta format).
  --hhself            Provide an existing HHsuite database of clusters to use.
  --hhdata            The path to the HHsuite data folder. If not provided will fetch
                      from the $HHLIB environment variable.
  --hhpfam            The pfam database formatted as an HHsuite database.
  --hhscop            The scop database formatted as an HHsuite database.
  --hhpdb             The pdb database formatted as an HHsuite database.
  --hhuniref          A uniref (or similar) database formatted as an HHsuite database.
                      Uniclust is a good preformatted option.
  --hhmatches_self    Precomputed matches of cluster HMMs against themselves, as a
                      directory containing an ffindex database of HHR results.
                      Files should be named as `db_hhr.ff{data,index}`.
  --hhmatches_pfam    Precomputed matches of cluster HMMs against Pfam.
                      Format is as for `--hhmatches_self`
  --hhmatches_scop    Precomputed matches of cluster HMMs against SCOP.
                      Format is as for `--hhmatches_self`
  --hhmatches_pdb     Precomputed matches of cluster HMMs against PDB.
                      Format is as for `--hhmatches_self`
  --hhmatches_uniref  Precomputed matches of cluster HMMs against uniref/uniclust.
                      Format is as for `--hhmatches_self`
```
