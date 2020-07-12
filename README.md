# Tools for working with mutation-specific primers
### used for 
+ validation of sequencing projects
+ running targeted deep-sequencing runs

### available as code best used from jupyter notebooks or docker container
+ #### for jupyter, add code to pythonpath and import tools (see next section):   
   * `from p3 import <tool>`
+ #### for docker:
   * `docker run -p <local-IP>:8888 -v $(pwd):/home/martin/work -v <your_static_folder>:/home/martin/static martin37szyska/primertools:<TAG>`
   * store your own code in `<mounted_volume>/code` and import directly (folder code automatically added to PYTHONPATH)
   * example notebooks are available in nb_templates

## The Tool kit
### run_primer3 - Retrieve primer sequences for given mutations
+ run_primer3 requires a pandas df containing mutations in the colunns Chr, Start, End, Ref, Alt
+ returns df populated with fwd/rev-Primer pairs and accessory info like insert size, annealing temp and insert sequence with mutation included


### Based on Jupyter Docker Stacks with modifications focused on size and simplicity

visit their documentation for more great content
* [Jupyter Docker Stacks on ReadTheDocs](http://jupyter-docker-stacks.readthedocs.io/en/latest/index.html)

Alpine versions are also available based on [alpine-miniconda3](https://hub.docker.com/r/frolvlad/alpine-miniconda3) with even smaller size

+ build the alpine version with:
    * `--build-arg=latest-alpine --build-arg=CURRENTUSER=10151`
