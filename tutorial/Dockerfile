# Dockerfile for building Atmospheric-Physics-Notebook container

FROM jupyter/minimal-notebook

LABEL maintainer="Markus Petters <mdpetter@ncsu.edu>"

USER root

# Install system packages dependencies
RUN apt-get update && \
    apt-get install -y --no-install-recommends \
    libfftw3-dev \
    vim \
    ffmpeg && \
    rm -rf /var/lib/apt/lists/*

# Environment Variables
ENV JULIA_DEPOT_PATH=/opt/julia
ENV JULIA_PKGDIR=/opt/julia
ENV JULIA_VERSION=1.4.2
ENV JULIA_PROJECT=$HOME

# Download and install julia version
RUN mkdir /opt/julia-${JULIA_VERSION} && \
    cd /tmp && \
    wget -q https://julialang-s3.julialang.org/bin/linux/x64/`echo ${JULIA_VERSION} | cut -d. -f 1,2`/julia-${JULIA_VERSION}-linux-x86_64.tar.gz && \
    tar xzf julia-${JULIA_VERSION}-linux-x86_64.tar.gz -C /opt/julia-${JULIA_VERSION} --strip-components=1 && \
    rm /tmp/julia-${JULIA_VERSION}-linux-x86_64.tar.gz
RUN ln -fs /opt/julia-*/bin/julia /usr/local/bin/julia

RUN mkdir /etc/julia && \
    echo "push!(Libdl.DL_LOAD_PATH, \"$CONDA_DIR/lib\")" >> /etc/julia/juliarc.jl && \
    mkdir $JULIA_PKGDIR && \
    chown $NB_USER $JULIA_PKGDIR && \
    chown $NB_USER  /opt/julia-${JULIA_VERSION} && \
    fix-permissions $JULIA_PKGDIR /opt/julia-${JULIA_VERSION}

USER $NB_UID

 # Download notebooks
RUN git clone  https://github.com/mdpetters/Data-Inversion-Tutorial.git && \
        cp -r $HOME/Data-Inversion-Tutorial/*.* . && \
        cp -r $HOME/Data-Inversion-Tutorial/* . &&  \
        rm -rf $HOME/work && \
        rm -rf $HOME/Data-Inversion-Tutorial

# Activate julia environment and precompile
RUN julia -e 'using Pkg; Pkg.instantiate()' && \
    julia -e 'using Pkg; Pkg.status()' && \
    julia -e 'using IJulia; IJulia.installkernel("Julia", "--depwarn=no")' && \
    julia -e 'using Pkg; Pkg.precompile()' && \
    julia -e 'println(Base.active_project())' 

RUN mv $HOME/.local/share/jupyter/kernels/julia* $CONDA_DIR/share/jupyter/kernels/ && \
    chmod -R go+rx $CONDA_DIR/share/jupyter && \
    rm -rf $HOME/.local && \
    fix-permissions $JULIA_PKGDIR $CONDA_DIR/share/jupyter

#USER root

# Install WebIO jupyter extension
RUN julia -e 'using WebIO; WebIO.install_jupyter_nbextension();' 

# Copy libraries for Fezzik precompile to succeed
USER root

RUN	cp $HOME/bootstrap.jl $JULIA_PKGDIR/packages/Fezzik/SfTjP/src/ && \
	chmod a+w ${JULIA_DEPOT_PATH}-${JULIA_VERSION}/lib/julia/ && \
    chmod a+w ${JULIA_DEPOT_PATH}-${JULIA_VERSION}/etc/julia/startup.jl 


USER $NB_UID

RUN echo 'using Fezzik; Fezzik.trace();' >> ${JULIA_DEPOT_PATH}-${JULIA_VERSION}/etc/julia/startup.jl && \
    jupyter nbconvert --to notebook --execute --ExecutePreprocessor.timeout=600 "Session 1 - Introduction.ipynb" --stdout >/dev/null && \ 
    jupyter nbconvert --to notebook --execute --ExecutePreprocessor.timeout=600 "Session 2 - The Forward Model.ipynb" --stdout >/dev/null && \ 
    jupyter nbconvert --to notebook --execute --ExecutePreprocessor.timeout=600 "Session 3 - Tikhonov Regularization.ipynb" --stdout >/dev/null && \ 
    jupyter nbconvert --to notebook --execute --ExecutePreprocessor.timeout=600 "Session 4 - Hands on Examples.ipynb" --stdout >/dev/null && \ 
    jupyter nbconvert --to notebook --execute --ExecutePreprocessor.timeout=600 "Session 5 - Summary and Perspective.ipynb" --stdout >/dev/null 
    
RUN julia -e 'using Fezzik; Fezzik.brute_build_julia(;clear_traces = true);'  
