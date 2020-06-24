################## BASE IMAGE ######################
FROM rocker/r-ubuntu:18.04   

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
