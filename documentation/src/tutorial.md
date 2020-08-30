# Tutorial

A tutorial using this package is available. The tutorial uses interactive jupyter notebooks,including audio and video clips. It features an introduction to differential mobility analyzers and the data inversion techniques. The tutorial is made available through Docker/DockerHub. This ensures low latency, no clutter from warnings, and a simple install process. 

## Docker
Install the [docker engine](https://docs.docker.com/get-docker). Then run the tutorial 
using the following command:

```bash
docker run -it -p 8888:8888 mdpetters/data-inversion-tutorial:v2009
```

This will download the container (~3GB) and exectute the container. It only needs to be downloaded once. Running the command will produce an http link similar to this one:

```bash
http://127.0.0.1:8888/?token=93b5e33a61654afded5f692325ac43f5c73eb6c58435196f
```

Note that the token is unique to each instance of the container. Follow the link and open the notebook:

```
Session 1 - Introduction.ipynb
```

Further instructions are printed inside the notebook. The tutorial takes ~1.5 hr to complete.

## Tutorial Files

The tutorial files are hosted in a separate [repository](https://github.com/mdpetters/Data-Inversion-Tutorial).