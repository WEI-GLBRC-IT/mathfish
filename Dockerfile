# mathFISH Docker Container
# Combines PHP/Apache, Python 3, ClustalW2, and UNAFold
# Target: linux/amd64

# PHP 8.4 on Debian bookworm. The minor-version tag picks up 8.4 patch
# releases on rebuild (CI builds with --pull); the floor check below fails
# the build on anything older than the PHP2Shell fix.
# Bookworm, not trixie: its gcc 12 is more tolerant of UNAFold 3.3's old C.
FROM php:8.4-apache-bookworm

LABEL maintainer="WEI IS Team"
LABEL description="mathFISH - FISH probe thermodynamic calculator"

# PHP2Shell (CVE-2026-109618, CVE-2026-109619) is fixed in 8.4.27.
RUN php -r 'exit(version_compare(PHP_VERSION, "8.4.27", "<") ? 1 : 0);' \
    || { echo "PHP $(php -r 'echo PHP_VERSION;') predates the PHP2Shell fix (8.4.27)"; exit 1; }

# Ensure /usr/local/bin is in PATH for all processes (where clustalw2 and UNAFold install)
ENV PATH="/usr/local/bin:/usr/local/sbin:${PATH}"

# Install system dependencies.
# numpy and matplotlib come from Debian: bookworm's pip refuses system-wide
# installs (PEP 668), and the old pip pins have no Python 3.11 builds.
RUN apt-get update && apt-get install -y --no-install-recommends \
    python3 \
    python3-numpy \
    python3-matplotlib \
    build-essential \
    perl \
    curl \
    && rm -rf /var/lib/apt/lists/*

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

# PHP: production defaults (display_errors Off, log_errors On) plus mathFISH overrides
RUN cp "$PHP_INI_DIR/php.ini-production" "$PHP_INI_DIR/php.ini"
COPY docker/php/zz-mathfish.ini "$PHP_INI_DIR/conf.d/"

# Apache: no directory listings, no version banners, no .htaccess
COPY docker/apache/zz-mathfish.conf /etc/apache2/conf-available/
RUN a2enconf zz-mathfish

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
