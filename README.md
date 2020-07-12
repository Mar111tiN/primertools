# Retrieve primer sequences for given mutations

+ python code is ready with `from p3 import run_primer3`

+ see example notebook in nb_templates

+ store your own code in `<mounted_volume>/code` and import directly (is added to PYTHONPATH)

+ you have to mount an additional folder containing chromosome-split genome-fasta-files

+ run with `docker run -p <local-IP>:8888 -v $(pwd):/home/martin/work -v <your_static_folder>:/home/martin/static martin37szyska/primertools:<TAG>`



### Based on Jupyter Docker Stacks with modifications focused on size and simplicity

visit their documentation for more great content
* [Jupyter Docker Stacks on ReadTheDocs](http://jupyter-docker-stacks.readthedocs.io/en/latest/index.html)

Alpine versions are also available based on [alpine-miniconda3](https://hub.docker.com/r/frolvlad/alpine-miniconda3) with even smaller size

+ build the alpine version with:
    * `--build-arg=latest-alpine --build-arg=CURRENTUSER=10151`