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
    libblas-dev \
    liblapack-dev \
    libatlas-base-dev \
    gfortran \
    gridengine-client \
    gridengine-common \
    gridengine-exec \
    make  \
    nano \
    python-matplotlib \
    python-skimage \
    supervisor \
    && \
    make -C /usr/local/src/gdal-docker install clean && \
    apt-get purge -y make

RUN curl https://bootstrap.pypa.io/get-pip.py -O
RUN python get-pip.py

COPY requirements.txt /tmp/

RUN pip install -r /tmp/requirements.txt

# Execute the antares container as root
USER root

# Output version and capabilities by default.
CMD gdalinfo --version && gdalinfo --formats && ogrinfo --formats
