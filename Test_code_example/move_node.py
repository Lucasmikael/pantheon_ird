import numpy as np
import matplotlib
import matplotlib as m
import matplotlib.path as mpath
import matplotlib.patches as mpatches
import matplotlib.pyplot as plt
import matplotlib.mlab as mlab
try:
    import matplotlib.pyplot as plt
except:
    raise

import networkx as nx
import math
matplotlib.use("TkAgg")

# try:
#     from networkx import graphviz_layout
#     layout=nx.graphviz_layout
# except ImportError:
#     print ("PyGraphviz not found; drawing with spring layout; will be slow.")
#     layout=nx.spring_layout


n = 50
p = 0.15
G = nx.MultiDiGraph()
G.add_node("plop")
G.add_node("plup")
G.add_node("plip")
G.add_edge("plop","plup")
G.add_edge("plop","plip")
#Gcc=nx.connected_component_subgraphs(G2)
# print (G.edges())
pos=nx.spring_layout(G, iterations=750)

Path = mpath.Path

fig = plt.figure()
ax = plt.gca()

pathdata1 = {}
patch1 = {}
n=-1
for u,v in G.edges():
      print (pos[u], pos[v])
      n = n +1
      pathdata1[n] = [
     (Path.MOVETO, pos[u]),    (Path.MOVETO, pos[v]),
      ]

      codes1, verts1 = zip(*pathdata1[n])
      path1 = mpath.Path(verts1, codes1)
      patch1[n] = mpatches.PathPatch(path1, facecolor='w', edgecolor='w', alpha=0.0)

      ax.add_patch(patch1[n])


class PathInteractor:
    """
    An path editor.

    Key-bindings

      't' toggle vertex markers on and off.  When vertex markers are on,
          you can move them, delete them


    """

    showverts = True
    epsilon = 5  # max pixel distance to count as a vertex hit

    def __init__(self, pathpatch):

        self.ax = pathpatch.axes
        canvas = self.ax.figure.canvas
        self.pathpatch = pathpatch
        self.pathpatch.set_animated(True)

        x, y = zip(*self.pathpatch.get_path().vertices)

        self.line, = ax.plot(x,y,color='gray',marker='o', markerfacecolor='r',animated=True)

        self._ind = None # the active vert

        canvas.mpl_connect('draw_event', self.draw_callback)
        canvas.mpl_connect('button_press_event', self.button_press_callback)
        canvas.mpl_connect('key_press_event', self.key_press_callback)
        canvas.mpl_connect('button_release_event', self.button_release_callback)
        canvas.mpl_connect('motion_notify_event', self.motion_notify_callback)
        self.canvas = canvas


    def draw_callback(self, event):
        self.background = self.canvas.copy_from_bbox(self.ax.bbox)
        self.ax.draw_artist(self.pathpatch)
        self.ax.draw_artist(self.line)
        self.canvas.blit(self.ax.bbox)

    def pathpatch_changed(self, pathpatch):
        'this method is called whenever the pathpatchgon object is called'
        # only copy the artist props to the line (except visibility)
        vis = self.line.get_visible()
        Artist.update_from(self.line, pathpatch)
        self.line.set_visible(vis)  # don't use the pathpatch visibility state


    def get_ind_under_point(self, event):
        'get the index of the vertex under point if within epsilon tolerance'

        # display coords
        xy = np.asarray(self.pathpatch.get_path().vertices)
        xyt = self.pathpatch.get_transform().transform(xy)
        xt, yt = xyt[:, 0], xyt[:, 1]
        d = np.sqrt((xt-event.x)**2 + (yt-event.y)**2)
        ind = d.argmin()

        if d[ind]>=self.epsilon:
            ind = None

        return ind

    def button_press_callback(self, event):
        'whenever a mouse button is pressed'
        if not self.showverts: return
        if event.inaxes==None: return
        if event.button != 1: return
        self._ind = self.get_ind_under_point(event)

    def button_release_callback(self, event):
        'whenever a mouse button is released'
        if not self.showverts: return
        if event.button != 1: return
        self._ind = None

    def key_press_callback(self, event):
        'whenever a key is pressed'
        if not event.inaxes: return
        if event.key=='t':
            self.showverts = not self.showverts
            self.line.set_visible(self.showverts)
            if not self.showverts: self._ind = None

        self.canvas.draw()

    def motion_notify_callback(self, event):
        'on mouse movement'
        if not self.showverts: return
        if self._ind is None: return
        if event.inaxes is None: return
        if event.button != 1: return
        x,y = event.xdata, event.ydata

        vertices = self.pathpatch.get_path().vertices

        vertices[self._ind] = x,y
        self.line.set_data(zip(*vertices))

        self.canvas.restore_region(self.background)
        self.ax.draw_artist(self.pathpatch)
        self.ax.draw_artist(self.line)
        self.canvas.blit(self.ax.bbox)


for n in range(len(G.edges())):
 interactor = PathInteractor(patch1[n])



plt.axis('off')
plt.show()
