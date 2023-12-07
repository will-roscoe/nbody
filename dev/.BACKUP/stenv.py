class staticPlane:
    def __init__(self,const_axis='z',axis_val=0.,color=('xkcd:azure',0.5)) -> None:
        self.c,self.v,self.col = *typecheck(((const_axis, str),(axis_val, NumType))),color
    
    
    def draw(self,ax,lims):
        pl_const,(xl,yl,zl) = np.array([[self.v,self.v],[self.v,self.v]]),lims
        
        points = {'x':(pl_const,np.array([[yl[0],yl[1]],[yl[0],yl[1]]]),
                        np.array([[zl[0],zl[0]],[zl[1],zl[1]]])),
                    'y':(np.array([[xl[0],xl[1]],[xl[0],xl[1]]]),pl_const,
                        np.array([[zl[0],zl[0]],[zl[1],zl[1]]])),
                    'z':(np.array([[xl[0],xl[1]],[xl[0],xl[1]]]),
                        np.array([[yl[0],yl[0]],[yl[1],yl[1]]]),pl_const)}
        ax.plot_surface(*points[self.c],zorder=1,color=self.col,clip_on=False) 


class staticVector:
    def __init__(self,vector):
        self.vec = _V(vector)
    
    
    def draw(self,ax):
        return


class staticEnv:
    def __init__(self,init_objs=list()):
        print(init_objs)
        self.objs = init_objs
    
    
    def add(self,object):
        typecheck((object,(*StaticObject,*Iterable)))
        if isinstance(object,staticEnv):
            eval(self.objs.append(ob) for ob in object.objs)
        
        elif isinstance(object,Iterable):
            eval(self.objs.append(ob) for ob in object if isinstance(ob,StaticObject))
        else:
            self.objs.append(object)
    
    
    def staticCube(self,dist=1.,color=('xkcd:azure',0.5)):
        typecheck((dist, NumType))
        self.add(list((staticPlane(ax,num,color) for num in (dist,-dist)) for ax in ('x','y','z')))
    
    
    def __getitem__(self, ind):
        _get_lookup = {**dict.fromkeys(['planes',staticPlane],list(plane for plane in self.objs if isinstance(plane,staticPlane))), 
                       **dict.fromkeys(['vectors',staticVector],list(vect for vect in self.objs if isinstance(vect,staticVector)))}
        return _get_lookup[ind]


    def draw(self, axes_obj):
        ax = typecheck((axes_obj, Axes3D))
        lims = (ax.get_xlim(),ax.get_ylim(),ax.get_zlim())
        (obj.draw(ax,lims) for obj in self.objs)
        
StaticObject =(staticPlane,staticVector,staticEnv)