FROM condaforge/miniforge3
WORKDIR /app/local_run
COPY local_run /app/local_run
RUN conda env create -f environment.yml && conda clean -afy
RUN conda install -n ppgwas_test -y jupyterlab && conda clean -afy
ENV PATH=/opt/conda/envs/ppgwas_test/bin:/opt/conda/bin:$PATH
RUN chmod +x ppgwas.sh
EXPOSE 8888
CMD bash -lc "jupyter lab --ip=0.0.0.0 --no-browser --NotebookApp.token='' --NotebookApp.password=''"
