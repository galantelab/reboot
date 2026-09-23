FROM rocker/r-ver:4.4.0

RUN apt-get update \
    && apt-get install -y --no-install-recommends \
      python3 \
    && rm -rf /var/lib/apt/lists/*

WORKDIR /app
COPY . /app

RUN R -e "install.packages('remotes')" \
    && R -e "remotes::install_local('.', upgrade = 'never', dependencies = TRUE)"

RUN chmod +x /app/inst/scripts/Reboot.R \
    && ln -s /app/inst/scripts/Reboot.R /usr/local/bin/Reboot.R

CMD ["Reboot.R"]
