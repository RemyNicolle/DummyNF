#non local based
FROM debian:10-slim
# FROM ubuntu:20.04
# FROM alpine:3.12
LABEL maintainer="Remy Nicolle (remy.c.nicolle@gmail.com)"


ENV STAR_VERSION=2.7.8a
ENV SUBREAD_VERSION=2.0.1
# ENV R_VERSION=4.0.4
ENV KALLISTO_VERSION=0.46.1
ENV SALMON_VERSION=1.4.0








##Requirements
# RUN yum -y update && yum install -y \
# curl-devel \
# ; yum clean all


# RUN yum -y update && yum install -y \
# gcc-c++ \
# gcc-gfortran  ; yum clean all \
# &&  cd /tmp \
# && curl -O https://cloud.r-project.org/src/base/R-4/R-${R_VERSION}.tar.gz \
# && tar -zxf R-${R_VERSION}.tar.gz \
# && cd /tmp/R-${R_VERSION} \
# && ./configure --with-readline=no --with-x=no \
# && make && make install \




# RUN  apk update && apk add --no-cache  R R-dev wget bash
# boost-dev xz-dev cmake curl-dev curl build-base




RUN  apt-get -y update \
&& apt-get -y  install r-base r-base-dev wget libboost-all-dev procps  \
liblzma-dev cmake libcurl4-openssl-dev curl \
&& apt-get autoremove && apt-get autoclean





RUN wget -nv https://github.com/COMBINE-lab/salmon/archive/v${SALMON_VERSION}.tar.gz \
&& tar -zxf v${SALMON_VERSION}.tar.gz \
&& cd salmon-${SALMON_VERSION}  \
&& mkdir build && cd build \
&& cmake .. -DCMAKE_INSTALL_PREFIX=/usr/local && make && make install




RUN cd /tmp \
&& wget  -nv https://github.com/alexdobin/STAR/archive/${STAR_VERSION}.tar.gz \
&& tar -zxf  ${STAR_VERSION}.tar.gz \
&& cp STAR-${STAR_VERSION}/bin/Linux_x86_64_static/STAR /usr/local/bin \
&& wget  -nv -O subread.tar.gz https://sourceforge.net/projects/subread/files/subread-${SUBREAD_VERSION}/subread-${SUBREAD_VERSION}-Linux-x86_64.tar.gz/download\
&& tar -zxf subread.tar.gz \
&& cp subread-${SUBREAD_VERSION}-Linux-x86_64/bin/featureCounts /usr/local/bin \
&& wget -nv https://github.com/pachterlab/kallisto/releases/download/v${KALLISTO_VERSION}/kallisto_linux-v${KALLISTO_VERSION}.tar.gz \
&& tar -zxf kallisto_linux-v${KALLISTO_VERSION}.tar.gz \
&& cp kallisto/kallisto /usr/local/bin \
&& chmod +x /usr/local/bin/* \
&& rm -rf /tmp/* ;

# COPY scripts /tmp/scripts
# HERE COPY SCRIPTS

# RUN cp /tmp/scripts/indexprepSTAR.sh /usr/local/bin/indexSTAR \
# && cp /tmp/scripts/countLexoFWD.sh  /usr/local/bin/countLexoFWD \
# RUN chmod +x /usr/local/bin/* \
# && rm -rf /tmp/* ;





# docker build  -t rnaseq:0.1 .
# docker run rnaseq:0.1  salmon -h

# docker tag rnaseq:0.1 gebican/rnaseq:0.1
# docker push gebican/rnaseq:0.1




# docker build --no-cache -t rnaseq:0.1 .

# docker container prune
# docker image prune



# Test index only
# docker run rnaseq:0.1  indexSTAR ftp.ensembl.org/pub/release-101/fasta/homo_sapiens/dna/Homo_sapiens.GRCh38.dna.primary_assembly.fa.gz ftp.ensembl.org/pub/release-101/gtf/homo_sapiens/Homo_sapiens.GRCh38.101.gtf.gz ~/Download/tmp 4 BS


# docker run rnaseq:0.1  indexSTAR ftp.ensembl.org/pub/release-101/gtf/caenorhabditis_elegans/Caenorhabditis_elegans.WBcel235.101.gtf.gz ftp.ensembl.org/pub/release-101/fasta/caenorhabditis_elegans/dna/Caenorhabditis_elegans.WBcel235.dna.toplevel.fa.gz  ~/Download/tmp 4 BS


# && curl --output /tmp/starcountshell.sh https://raw.githubusercontent.com/gebican/RNASeq_container/main/QuantSeqSE/STARnCount.sh?token=ABAHQVZWGTZRFEDH2FML5P27YT4L6 \
# && cp /tmp/starcountshell.sh /usr/local/bin/MAPnCOUNT \
