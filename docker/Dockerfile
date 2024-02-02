FROM mambaorg/micromamba

# for building the conda env
COPY --chown=$MAMBA_USER:$MAMBA_USER environment.yml /tmp/environment.yml
# Build the environment based on the yml file
RUN micromamba install -y -n base -f /tmp/environment.yml \
  && micromamba clean --all --yes

ADD create_gene_mappings.R /usr/local/bin/
ADD map.R /usr/local/bin/
ADD run.sh /usr/local/bin/

USER root
RUN mkdir -p /opt/resources && \
  chmod a+r /opt/resources && \
  chown $MAMBA_USER:$MAMBA_USER /opt/resources

# Run the install scripts and create the mapping files.
USER $MAMBA_USER
ARG MAMBA_DOCKERFILE_ACTIVATE=1
RUN Rscript /usr/local/bin/create_gene_mappings.R human /opt/resources/human_mappings.tsv
RUN Rscript /usr/local/bin/create_gene_mappings.R mouse /opt/resources/mouse_mappings.tsv

