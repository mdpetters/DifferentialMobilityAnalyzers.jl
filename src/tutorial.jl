function launch_tutorial()
    path = pathof(DifferentialMobilityAnalyzers)
    WebIO.install_jupyter_nbextension()
    
    dd = dirname(path)
    cd(dd)
    cd("../tutorial")
    notebook(dir = pwd())
    
end