import pandas

def outdeg_calculator(tfnet):
    return(tfnet.groupby(by="tf").sum())

def indeg_calculator(targetnet):
    return(targetnet.groupby(by="gene").sum())
