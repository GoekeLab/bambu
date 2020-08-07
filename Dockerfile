################## BASE IMAGE ######################
FROM rocker/r-ubuntu:20.04


ENV DEBIAN_FRONTEND noninteractive

RUN apt-get update && apt-get -y install --no-install-recommends --no-install-suggests \
        ca-certificates software-properties-common gnupg2 gnupg1 \
      && apt-key adv --keyserver keyserver.ubuntu.com --recv-keys E298A3A825C0D65DFD57CBB651716619E084DAB9 \
      && add-apt-repository 'deb https://cloud.r-project.org/bin/linux/ubuntu focal-cran40/' \
      && apt-get install r-base -y 

# Install required libraries -- using prebuild binaries where available
RUN apt-get update && apt-get install -y \
aptitude \
libcurl4-openssl-dev \
libxml2-dev \
git \
r-cran-devtools \
r-cran-git2r \
r-cran-xml \
r-cran-rcurl \
sudo

# Install plr -- for now (?) from GH; also install visualTest
RUN installGithub.r Goekelab/bambu \
&& rm -rf /tmp/downloaded_packages/
