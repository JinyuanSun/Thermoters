# Docker Usage

This repository can be used inside a Docker container. The provided `Dockerfile` installs all Python dependencies listed in `requirements.txt` and copies the repository files into the image.

## Build the image

```bash
docker build -t thermoters .
```

## Run the container

To start an interactive Python session:

```bash
docker run -it thermoters
```

You can also mount the repository directory to keep changes outside the container:

```bash
docker run -it -v $(pwd):/app thermoters
```

