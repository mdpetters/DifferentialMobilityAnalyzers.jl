# Tutorial

A tutorial using this package is available. The tutorial uses interactive jupyter notebooks,
including audio and video clips. It features an introduction to differential mobility
analyzers and the data inversion techniques. Users new to Julia should run this tutorial through 
Docker/DockerHub. This ensures low latency and a simple install process. 

## Docker
Install the docker engine from https://docs.docker.com/get-docker/. Then run the tutorial 
using the following command

```bash
docker run -it -p 8888:8888 mdpetters/data-inversion-tutorial:v2009
```

Running the command will produce an http link similar to this one:

```bash
http://127.0.0.1:8888/?token=93b5e33a61654afded5f692325ac43f5c73eb6c58435196f
```

Note that the token is unique to each instance of the container. Follow the link and open 
the notebook:

```
Session 1 - Introduction.ipynb
```

Further instructions are printed inside the notebook. The tutorial takes ~1.5 hr to complete.

## Local

The notebooks are also included locally in the tutorial folder and can be executed through 
IJulia.

