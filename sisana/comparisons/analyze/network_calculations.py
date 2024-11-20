import pandas

def outdeg_calculator(tfnet):
    return(tfnet.groupby(by="TF").sum())

def indeg_calculator(targetnet):
    return(targetnet.groupby(by="Target").sum())
