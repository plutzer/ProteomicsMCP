# Docker Deployment for CPTAC Proteomics Explorer

This guide explains how to build and run the CPTAC Proteomics Explorer using Docker.

## Quick Start

### Using Docker Compose (Recommended)

```bash
# Build and start
docker-compose up -d

# View logs
docker-compose logs -f

# Stop the container
docker-compose down
```

### Using Docker CLI

```bash
# Build the image
docker build -t cptac-explorer .

# Run the container
docker run -d -p 3838:3838 --name cptac-explorer cptac-explorer

# View logs
docker logs -f cptac-explorer

# Stop
docker stop cptac-explorer && docker rm cptac-explorer
```

## Accessing the Application

Once the container is running, access the app at:
- **Local:** http://localhost:3838
- **Network:** http://YOUR_IP_ADDRESS:3838

## Container Details

### Base Image
- **Python 3.11-slim** - Lightweight Debian-based image

### Exposed Ports
- **3838** - Shiny web application

### Included Data
The container includes pre-downloaded CPTAC data for all 9 cancer types:
- BRCA (Breast Cancer)
- COAD (Colon Adenocarcinoma)
- CCRCC (Clear Cell Renal Cell Carcinoma)
- GBM (Glioblastoma)
- HNSCC (Head and Neck Squamous Cell Carcinoma)
- LSCC (Lung Squamous Cell Carcinoma)
- LUAD (Lung Adenocarcinoma)
- OV (Ovarian Cancer)
- PDAC (Pancreatic Ductal Adenocarcinoma)

## Troubleshooting

### Container Won't Start
Check logs for errors:
```bash
docker-compose logs
```

Common issues:
- Port 3838 already in use: Change port in `docker-compose.yml`
- Build failures: Ensure all files are present

### Slow Startup
The application needs to preprocess data at startup, which takes 5-10 minutes.
The healthcheck has a `start_period: 120s` to account for this.

## Development Workflow

### Rebuild After Code Changes

```bash
# Using Docker Compose
docker-compose down
docker-compose up --build -d

# Using Docker CLI
docker stop cptac-explorer && docker rm cptac-explorer
docker build -t cptac-explorer .
docker run -d -p 3838:3838 --name cptac-explorer cptac-explorer
```

### Interactive Shell Access

```bash
docker-compose exec cptac-explorer /bin/bash
```

## Production Deployment

For production environments, consider:

1. **Use a reverse proxy** (nginx, traefik) for HTTPS
2. **Add authentication** using nginx basic auth or OAuth2 proxy
3. **Resource limits:**
   ```yaml
   deploy:
     resources:
       limits:
         cpus: '2'
         memory: 4G
   ```

## File Structure in Container

```
/app/
├── local_data_loader.py      # Local data loading module
├── cptac_proteomics.py       # CPTAC MCP tools
├── GUI/
│   ├── cptac_explorer_app.py # Main Shiny app
│   └── cptac_backend.py      # Data processing backend
├── data/                     # Pre-downloaded CPTAC data
│   ├── bcm-brca/
│   ├── bcm-coad/
│   └── ...
└── Datasets/                 # PSP data files
```

## Cleaning Up

```bash
# Stop and remove containers
docker-compose down

# Remove images
docker rmi cptac-explorer

# Remove all unused Docker resources
docker system prune -a
```
