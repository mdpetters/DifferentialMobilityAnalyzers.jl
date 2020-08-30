function launch_tutorial()
    path = pathof(DifferentialMobilityAnalyzers)
    dd = dirname(path)
    cd(dd)
    cd("../tutorial")
    notebook(dir = pwd())
end