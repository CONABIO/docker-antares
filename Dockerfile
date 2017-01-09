##
# geodata/gdal
#
# This creates an Ubuntu derived base image that installs the latest GDAL
# subversion checkout compiled with a broad range of drivers.  The build process
# is based on that defined in
# <https://github.com/OSGeo/gdal/blob/trunk/.travis.yml>
#

# Ubuntu 14.04 Trusty Tahyr
FROM ubuntu:trusty

MAINTAINER Amaury Gutierrez <amaury.gtz@gmail.com>

# Install the application.
ADD . /usr/local/src/gdal-docker/
RUN apt-get update -y && \
    apt-get install -y \
    curl \
    git \
    gfortran \
    gridengine-client \
    gridengine-common \
    gridengine-exec \
    htop \
    libblas-dev \
    libeigen3-dev \
    liblapack-dev \
    libatlas-base-dev \
    libtbb-dev \
    make \
    nano \
    python-matplotlib \
    python-skimage \
    supervisor \
    wget \
    && \
    make -C /usr/local/src/gdal-docker install clean && \
    make -C /usr/local/src/gdal-docker/correlation && \
    apt-get purge -y make
RUN ln -s /usr/local/src/gdal-docker/correlation/correlation /usr/local/bin/correlation
RUN curl https://bootstrap.pypa.io/get-pip.py -O
RUN python get-pip.py

COPY requirements.txt /tmp/

RUN pip install -r /tmp/requirements.txt

# Execute the antares container as root
USER root

RUN export MASTER_HOSTNAME=master
RUN echo "gridengine-common       shared/gridenginemaster string  $MASTER_HOSTNAME" | sudo debconf-set-selections
RUN echo "gridengine-common       shared/gridenginecell   string  default" | sudo debconf-set-selections
RUN echo "gridengine-common       shared/gridengineconfig boolean false" | sudo debconf-set-selections
RUN echo "gridengine-client       shared/gridenginemaster string  $MASTER_HOSTNAME" | sudo debconf-set-selections
RUN echo "gridengine-client       shared/gridenginecell   string  default" | sudo debconf-set-selections
RUN echo "gridengine-client       shared/gridengineconfig boolean false" | sudo debconf-set-selections
RUN echo "postfix postfix/main_mailer_type        select  No configuration" | sudo debconf-set-selections
RUN service postfix stop
RUN update-rc.d postfix disable
RUN service gridengine-exec restart

# Output version and capabilities by default.
CMD gdalinfo --version && gdalinfo --formats && ogrinfo --formats
