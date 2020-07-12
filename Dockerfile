ARG TAG=3.8

FROM martin37szyska/jupyter-minimal:$TAG

ARG CURRENTUSER=$NB_USER

USER root

# install primer3
RUN pip install --upgrade pip && \
    pip install --no-cache-dir primer3-py && \
    fix-permissions /home/$NB_USER

# add code_master to the python path and a possibly local folder code to PYTHONPATH
RUN echo "sys.path.insert(0,'/home/martin/code_master')" >> .ipython/profile_default/startup/startup.py && \
    echo "sys.path.insert(0,'/home/martin/work/code')" >> .ipython/profile_default/startup/startup.py

# copy the primer3-manual for reference
COPY manual ${HOME}/manual

# copy the registry code to code_master
COPY code ${HOME}/code_master

# copy the notebook templates to nb_templates for local use
COPY nb ${HOME}/nb_templates

COPY testdata ${HOME}/testdata


# Switch back to martin to avoid accidental container runs as root
USER $CURRENTUSER