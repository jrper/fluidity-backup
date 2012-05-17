import vtktools

def contour2D(fname,vname,values,**kwargs):

    import pylab as p

    file=vtktools.vtu(fname)

    pts=file.GetLocations()

    x=p.linspace(min(pts[:,0]),max(pts[:,0]),5001)
    y=p.linspace(min(pts[:,1]),max(pts[:,1]),5001)

               
    X,Y=p.meshgrid(x,y)

    grd=p.array((X.ravel(),Y.ravel(),0*X.ravel())).T

    Z=file.ProbeData(grd,vname).reshape(X.shape)

    p.figure()
    p.contour(X,Y,Z,values)


    return X,Y,Z
