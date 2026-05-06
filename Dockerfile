# mathFISH Docker Container
# Combines PHP/Apache, Python 3, ClustalW2, and UNAFold
# Target: linux/amd64

FROM php:7.4-apache

LABEL maintainer="WEI IS Team"
LABEL description="mathFISH - FISH probe thermodynamic calculator"

# Ensure /usr/local/bin is in PATH for all processes (where clustalw2 and UNAFold install)
ENV PATH="/usr/local/bin:/usr/local/sbin:${PATH}"

# Install system dependencies
RUN apt-get update && apt-get install -y \
    python3 \
    python3-pip \
    python3-dev \
    build-essential \
    gfortran \
    perl \
    curl \
    && rm -rf /var/lib/apt/lists/*

# Install Python packages (pinned versions for compatibility)
RUN pip3 install --no-cache-dir \
    matplotlib==3.3.4 \
    numpy==1.19.5 \
    Pillow==8.4.0 \
    python-dateutil==2.9.0.post0 \
    pyparsing==3.1.2

# Copy and build ClustalW2
COPY dependencies/clustalw-2.1.tar.gz /tmp/
RUN cd /tmp && \
    tar -xzf clustalw-2.1.tar.gz && \
    cd clustalw-2.1 && \
    ./configure && \
    make && \
    make install && \
    chmod a+rx /usr/local/bin/clustalw2 && \
    cd / && rm -rf /tmp/clustalw-2.1*

# Copy and build UNAFold
COPY dependencies/unafold-3.3.tar.gz /tmp/
RUN cd /tmp && \
    tar -xzf unafold-3.3.tar.gz && \
    cd unafold-3.3 && \
    ./configure && \
    make && \
    make install && \
    cd / && rm -rf /tmp/unafold-3.3*

# Configure Apache - set document root and enable modules
RUN a2enmod rewrite
ENV APACHE_DOCUMENT_ROOT=/var/www/html

# Set working directory (PHP scripts use relative paths to userFiles/)
WORKDIR /var/www/html

# Copy application files (inputFiles/ is inside web/)
COPY --chown=www-data:www-data web/ /var/www/html/

# Create userFiles directory with proper permissions
# www-data (uid 33) is the Apache user in Debian
RUN mkdir -p /var/www/html/userFiles && \
    chown -R www-data:www-data /var/www/html && \
    chmod -R 755 /var/www/html && \
    chmod 2775 /var/www/html/userFiles

# Expose port 80
EXPOSE 80

# Health check
HEALTHCHECK --interval=30s --timeout=3s \
    CMD curl -f http://localhost/ || exit 1

# Start Apache in foreground
CMD ["apache2-foreground"]
