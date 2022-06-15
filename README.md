# Tools for working with mutation-specific primers
### used for 
+ validation of sequencing projects
+ running targeted deep-sequencing runs

### available as code best used from jupyter notebooks or docker container
+ #### for jupyter:
   * clone the repository into \<your folder\> and move to into primertools:
      + `cd <your_folder> && git clone https://github.com/Mar111tiN/primertools.git && cd primertools`
   * create conda environment to run the notebooks (for AppleSilicon use env/primer3_M1-env.yml):
      + `conda env create -n primer3-env -f env/primer3-env.yml`
   * run the jupyter notebook and try primer 3 on your testdata with primer3_test.ipynb
      + instructions are included in the notebook 

+ #### for docker:
   * see on [Docker Hub](https://hub.docker.com/repository/docker/martin37szyska/primertools)
   * `docker run -p <local-IP>:8888 -v $(pwd):/home/martin/work -v <your_static_folder>:/home/martin/static martin37szyska/primertools:<TAG>`
   * store your own code in `<mounted_volume>/code` and import directly (folder code automatically added to PYTHONPATH)
   * example notebooks are available in nb_templates

## The Tool kit
### run_primer3 - Retrieve primer sequences for given mutations
+ run_primer3 requires a pandas dataframe containing mutations in the columns Chr, Start, End, Ref, Alt
+ returns df populated with fwd/rev-Primer pairs and accessory info like insert size, annealing temp and insert sequence with mutation included

### check_primerDB - search for suitable primers in primer list
+ check_primerDB requires mutation dataframe of same format as run_primer3
+ returns tuple of 
   * same mutation dataframe with added DBhits column giving the number of hits in the DB
   * dataframe containing the primers from the database for each mutation with a hit
+ try out in notebook primer3-test.ipynb 

### Docker image is based on modified Jupyter Docker Stacks with modifications focused on size and simplicity

visit their documentation for more great content
* [Jupyter Docker Stacks on ReadTheDocs](http://jupyter-docker-stacks.readthedocs.io/en/latest/index.html)
