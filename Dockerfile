FROM condaforge/miniforge3
SHELL ["/bin/bash","-lc"]
WORKDIR /app/local_run
COPY local_run/environment.yml /tmp/environment.yml
RUN conda install -n base -c conda-forge -y mamba && mamba env create -f /tmp/environment.yml && conda clean -afy
ENV PATH=/opt/conda/envs/ppgwas_test/bin:/opt/conda/bin:$PATH
RUN conda install -n ppgwas_test -y jupyterlab && conda clean -afy
COPY local_run /app/local_run
RUN chmod +x ppgwas.sh
EXPOSE 8888
CMD jupyter lab --ip=0.0.0.0 --no-browser --NotebookApp.token='' --NotebookApp.password=''

