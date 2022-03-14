From r-base:4.1.2

RUN apt-get update &&\
  apt-get install -y libcurl4-openssl-dev \
      libxml2-dev \
      libssl-dev

# Install all software under /opt/software:
RUN mkdir -p /opt/software 
RUN mkdir -p /opt/software/resources 

# Now add the necessary scripts:
ADD install_packages.R /opt/software
ADD create_gene_mappings.R /opt/software/
ADD map.R /opt/software/
ADD run.sh /opt/software
RUN chmod +x /opt/software/run.sh

# Run the install scripts and create the mapping files.
RUN Rscript /opt/software/install_packages.R
RUN Rscript /opt/software/create_gene_mappings.R human /opt/software/resources/human_mappings.tsv
RUN Rscript /opt/software/create_gene_mappings.R mouse /opt/software/resources/mouse_mappings.tsv

ENTRYPOINT ["/bin/bash"]

