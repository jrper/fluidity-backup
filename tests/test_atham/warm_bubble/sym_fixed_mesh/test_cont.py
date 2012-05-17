import pylab as p
from operator import lt, gt

def cont_tri(pts,vals,v,**kwargs):
    if not 'color' in kwargs:
        kwargs['color']='b'
    if (sum(vals>v)==2 or sum(vals<v)==2):
        if lt(vals[0],v):
            op=gt
        else:
            op=lt
        t1=op(vals[1],v)
        t2=op(vals[2],v)
        if t1 and t2 :
            P1=pts[:][0]+(v-vals[0])*(pts[:][1]-pts[:][0])/(vals[1]-vals[0])
            P2=pts[:][0]+(v-vals[0])*(pts[:][2]-pts[:][0])/(vals[2]-vals[0])
        elif t1:
            P1=pts[:][0]+(v-vals[0])*(pts[:][1]-pts[:][0])/(vals[1]-vals[0])
            P2=pts[:][1]+(v-vals[1])*(pts[:][2]-pts[:][1])/(vals[2]-vals[1])
        else:
            P1=pts[:][0]+(v-vals[0])*(pts[:][2]-pts[:][0])/(vals[2]-vals[0])
            P2=pts[:][1]+(v-vals[1])*(pts[:][2]-pts[:][1])/(vals[2]-vals[1])
        p.plot([P1[0],P2[0]],[P1[1],P2[1]],**kwargs)

#pts=p.array([[0,0],[1,0],[0,1]])
#vals=p.array([0,1,1])
#cont_tri(pts,vals,0.5)


def superpositions(f):
    X=f.GetLocations()
    node_dict={}
    for k,r in enumerate(X):
        R=tuple(r)
        node_dict[R]=node_dict.get(R,[])+[k]
    return node_dict

def overlap(f,name,v):
    nd=superpositions(f)
    V=f.GetField(name)
    print V
    overlaps=p.zeros(len(V))

    for I in nd.values():
        print I
        flag=any(V[I]>v) and any(V[I]<v)
        for i in I:
            overlaps[i]=flag

    return overlaps    

def cont_vtk(f,name,v,**kwargs):
    X=f.GetLocations()
    V=f.GetField(name)

    overlaps=overlap(f,name,v)
    
    S=p.array([[0,3,5],[1,3,4],[2,4,5],[3,4,5]])
    l3=p.array([[0,1],[1,2],[2,0]])
    l6=p.array([[0,5],[5,2],[2,4],[4,1],[1,3],[3,0]])

    itr=0
    pts0=f.GetCellPoints(0)
    while True:
        pts=f.GetCellPoints(itr)
        if len(pts)>6:
            break
        if itr>0 and all(pts==pts0):
            break
        if len(pts)==3:
            cont_tri(X[pts],V[pts],v)
            for l in l3:
                if all(overlaps[pts[l]]):
                    p.plot([X[pts[l[0]]][0],X[pts[l[1]]][0]],
                           [X[pts[l[0]]][1],X[pts[l[1]]][1]])
            print itr
        else:
            for k,s in enumerate(S):
                cont_tri(X[pts[s]],V[pts[s]],v)
            print '%d:%d'%(itr,k)
            for l in l6:
                if all(overlaps[pts[l]]):
                    p.plot([X[pts[l[0]]][0],X[pts[l[1]]][0]],
                           [X[pts[l[0]]][1],X[pts[l[1]]][1]])
        itr+=1




    
    

        
