# mathFISH

Web application for calculating thermodynamic parameters of FISH (Fluorescence In Situ Hybridization) probes.

## Table of Contents

- [Design](#design)
  - [Components](#components)
  - [Why these versions are pinned](#why-these-versions-are-pinned)
  - [Repository layout](#repository-layout)
- [Install and Run](#install-and-run)
  - [Prerequisites](#prerequisites)
  - [Obtain the UNAFold tarball](#obtain-the-unafold-tarball)
  - [Build and run](#build-and-run)
  - [Updating an existing deployment](#updating-an-existing-deployment)
- [Configuration](#configuration)
- [Testing](#testing)
- [Troubleshooting](#troubleshooting)

## Design

mathFISH ships as a single Docker image. The image is built locally from a UNAFold tarball you supply — UNAFold's license forbids redistribution, so the tarball cannot be committed to this repository or distributed in a prebuilt image. Everything else (PHP, Apache, Python code, ClustalW, reference data) is in the image.

### Components

| Component  | Version | Source           | Purpose                                                  |
|------------|---------|------------------|----------------------------------------------------------|
| PHP        | 7.4     | `php:7.4-apache` | Web frontend                                             |
| Apache     | 2.4     | `php:7.4-apache` | Web server                                               |
| Python     | 3.9     | Debian bullseye  | Computation scripts                                      |
| numpy      | 1.19.5  | pip              | Numerical computing                                      |
| matplotlib | 3.3.4   | pip              | Plot generation                                          |
| Pillow     | 8.4.0   | pip              | Image processing                                         |
| ClustalW2  | 2.1     | `dependencies/`  | Pairwise sequence alignment (tarball in repo already)    |
| UNAFold    | 3.3     | user-supplied    | Nucleic acid folding (`hybrid-ss-min`, `hybrid-2s.pl`)   |

### Why these versions are pinned

The image builds from the `php:7.4-apache` upstream base, which fixes PHP to 7.4.33 and Apache to 2.4.54 (on Debian bullseye). PHP 7.4 is the last interpreter version mathFISH has been validated against — patterns used in the frontend (`preg_replace`, array indexing style, implicit `{$var}` interpolation) have been deprecated or removed in newer PHP, and the frontend has not been retested under 8.x. Python pins ride along: bullseye ships Python 3.9.2, and the `numpy 1.19.5` / `matplotlib 3.3.4` / `Pillow 8.4.0` triple is the pip-installable set that's ABI-compatible with that interpreter and produces the plots matching the upstream mathFISH publication. Bumping any of these should be done deliberately, with a full pass of [tests/TESTING.md](tests/TESTING.md) before trusting the result.

### Repository layout

```
mathfish/
├── Dockerfile               # builds the runtime image (PHP/Apache + Python + ClustalW + UNAFold)
├── docker-compose.yml       # one-service compose for the runtime container
├── readme.md                # this file
├── .dockerignore, .gitignore
├── dependencies/            # build inputs
│   ├── clustalw-2.1.tar.gz  # committed (permissive license, ships in Bioconda)
│   └── unafold-3.3.tar.gz   # NOT committed (licensed); placed by user before build
├── tests/                   # smoke tests + per-feature procedures
│   ├── TESTING.md
│   ├── run_cli_tests.sh
│   └── …test fixtures
└── web/                     # application code copied into the image at /var/www/html/
    ├── *.php                # frontend pages
    ├── *.py                 # computation scripts
    └── inputFiles/          # reference rRNA sequences and stacking parameters
```

## Install and Run

### Prerequisites

- Docker Engine (or Docker Desktop) with Docker Compose v2
- A copy of `unafold-3.3.tar.gz` you've obtained under your own academic license (see below)
- On Apple Silicon / other ARM hosts: Docker will emulate `linux/amd64`, slower but functional

### Obtain the UNAFold tarball

UNAFold's license does not permit redistribution, so this repository does not include it. Before building the image you must obtain a copy yourself:

1. Register at [unafold.org](https://www.unafold.org/) under your institution's academic license.
2. Download `unafold-3.3.tar.gz`.
3. Place it at `./dependencies/unafold-3.3.tar.gz`.

The Dockerfile compiles UNAFold from this tarball at build time and bakes the resulting binaries into the image. The tarball stays out of git via `.gitignore`, so it persists across `git pull` and won't be accidentally committed.

### Build and run

```bash
git clone <this repository> mathfish
cd mathfish

# place dependencies/unafold-3.3.tar.gz here (see above)

docker compose build
docker compose up -d

# open http://localhost
```

Defaults: container name `mathfish`, host port 80.

### Updating an existing deployment

```bash
cd /path/to/mathfish
git pull
docker compose build --no-cache && docker compose up -d
```

If the tarball was lost (for example the working tree was wiped), re-stage `dependencies/unafold-3.3.tar.gz` from your secure copy before building.

## Configuration

### Ports

Default: `80:80` (host:container). Edit [docker-compose.yml](docker-compose.yml) to change:

```yaml
ports:
  - "8080:80"   # host:container
```

### Environment variables

| Variable | Default | Description |
|---|---|---|
| `TZ` | `America/Chicago` | Container timezone |

### Volume mounts (optional)

For debugging, file persistence, or live editing, uncomment the volumes block in [docker-compose.yml](docker-compose.yml):

```yaml
volumes:
  - ./userFiles:/var/www/html/userFiles
```

Create the host directory first: `mkdir -p ./userFiles && chmod 777 ./userFiles`.

### Container management

```bash
# Start / stop
docker compose up -d
docker compose down

# Logs
docker logs mathfish

# Shell into container
docker exec -it mathfish bash

# Status / health
docker compose ps
```

## Testing

See [tests/TESTING.md](tests/TESTING.md) for:

- CLI tests runnable inside the Docker container
- Web UI test procedures for all 8 analysis features
- Expected output values for validation
- A copy-pasteable General Analysis smoke test

To run all 8 CLI cases against the running container and produce a single log:

```bash
docker cp tests/run_cli_tests.sh mathfish:/tmp/run_cli_tests.sh
docker exec -w /var/www/html mathfish bash /tmp/run_cli_tests.sh > tests.out
```

## Troubleshooting

Container won't start:
```bash
docker logs mathfish
```

Missing UNAFold tarball at build time — `docker compose build` errors with `COPY failed`:
```bash
ls dependencies/unafold-3.3.tar.gz   # must exist before building
```

Check binaries are on PATH:
```bash
docker exec mathfish which clustalw2 hybrid-ss-min python3
```

Check `userFiles/` permissions — should show `drwxrwsr-x` (2775) owned by `www-data:www-data`:
```bash
docker exec mathfish ls -la /var/www/html/userFiles/
```

Health check failing:
```bash
docker exec mathfish curl -f http://localhost/
```
