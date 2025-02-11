
#For non convex domain-------------------------------

def show(self):
    fig = go.Figure()
    self.addPlot(fig)
    fig.update_layout(yaxis=dict(
        scaleanchor='x',
        scaleratio=1))
    fig.show()

def addPlot(self, fig):
    fig.add_trace(go.Scatter(x=[P[0] for P in self.listOuterVertex]+[self.listOuterVertex[0][0]],
                            y=[P[1] for P in self.listOuterVertex]+[self.listOuterVertex[0][1]],
                            fill="toself", fillcolor="White", mode="lines"))
    for edge in self.wallHolesEdges:
        fig.add_shape(type="line",
                x0=self.wallHolesVertex[edge[0]][0],
                y0=self.wallHolesVertex[edge[0]][1],
                x1=self.wallHolesVertex[edge[1]][0],
                y1=self.wallHolesVertex[edge[1]][1],
                line=dict(
                    color="LightSeaGreen",
                    width=2,
                ))
    for exit in self.exitList:
        fig.add_shape(type="line",
                x0=exit[0][0],
                y0=exit[0][1],
                x1=exit[1][0],
                y1=exit[1][1],
                line=dict(
                    color="Red",
                    width=2,
                ))

#for mesh2D -------------------------------------

def show(self):
    fig = go.Figure()
    self.domain.addPlot(fig)
    for T in self.triangles:
        fig.add_trace(go.Scatter(x=[self.vertices[i][0] for i in T]+[self.vertices[T[0]][0]],
                                y=[self.vertices[i][1] for i in T]+[self.vertices[T[0]][1]],
                        fill="toself",
                        fillcolor="White",
                        mode="lines",
                        line=dict(
                            color="Black",
                            width=1
                         )))
    fig.update_layout(yaxis=dict(
        scaleanchor='x',
        scaleratio=1))
    fig.show()

def addPlot(self, fig):
    for T in self.triangles:
        fig.add_trace(go.Scatter(x=[self.vertices[i][0] for i in T]+[self.vertices[T[0]][0]],
                                y=[self.vertices[i][1] for i in T]+[self.vertices[T[0]][1]],
                        fill="toself",
                        fillcolor="White",
                        mode="lines",
                        line=dict(
                            color="Black",
                            width=1
                         )))

#for cell value map

def show(self):
    fig = go.Figure()
    #self.Mesh.domain.addPlot(fig)
    for j,T in enumerate(self.Mesh.triangles):
        fig.add_trace(go.Scatter(x=[self.Mesh.vertices[i][0] for i in T]+[self.Mesh.vertices[T[0]][0]],
                                y=[self.Mesh.vertices[i][1] for i in T]+[self.Mesh.vertices[T[0]][1]],
                        fill="toself",
                        hoverinfo = "none",
                        showlegend = False,
                        mode="none",
                        fillcolor ='rgb('+str( int(255*min(1,max(self.values[j],0))) )+',0,0)'
                        ))
    fig.update_layout(yaxis=dict(
        scaleanchor='x',
        scaleratio=1))
    fig.show()

def getScatter(self):
    L = []
    for j,T in enumerate(self.Mesh.triangles):
        L.append(go.Scatter(x=[self.Mesh.vertices[i][0] for i in T]+[self.Mesh.vertices[T[0]][0]],
                                y=[self.Mesh.vertices[i][1] for i in T]+[self.Mesh.vertices[T[0]][1]],
                        fill="toself",
                        hoverinfo = "none",
                        showlegend = False,
                        mode="none",
                        fillcolor ='rgb('+str( int(255*min(1,max(self.values[j],0))) )+',0,0)'
                        ))
    return L

# for vertex value map----------------------------

def show3D(self):
    fig = go.Figure()
    fig.add_trace(go.Mesh3d(x=[P[0] for P in self.Mesh.vertices],
                            y = [P[1] for P in self.Mesh.vertices],
                            z = self.values,
                            opacity=1,
                            color='rgba(244,22,100,0.6)'))
    fig.show()

def show(self, grid=False, colorscale_name = 'viridis'):
    fig = go.Figure()
    self.Mesh.domain.addPlot(fig)
    if(grid):
        self.Mesh.addPlot(fig)
    fig.add_trace(go.Scatter(x=[P[0] for P in self.Mesh.vertices],
                            y = [P[1] for P in self.Mesh.vertices],
                    hoverinfo = "none",
                    showlegend = False,
                    mode="markers",
                    marker = dict(
                    color = self.values,
                    colorscale = colorscale_name
                    )))
    fig.update_layout(yaxis=dict(
        scaleanchor='x',
        scaleratio=1))
    fig.show()


def computeGradientFlow(self):
    LTrianglesGrad = []
    for point in range(len(self.Mesh.vertices)):
        if(self.values[point] < 0):
            print("GROS PROBLEME");
    for triangle in self.Mesh.triangles: #On a des déterminants égaux à 0....
            det =  ((self.Mesh.vertices[triangle[1]][0] - self.Mesh.vertices[triangle[0]][0])*(self.Mesh.vertices[triangle[2]][1] - self.Mesh.vertices[triangle[0]][1])
                    - (self.Mesh.vertices[triangle[2]][0] - self.Mesh.vertices[triangle[0]][0])*(self.Mesh.vertices[triangle[1]][1] - self.Mesh.vertices[triangle[0]][1]))
            if(det == 0):
                print("det nul :")
                print(str(self.Mesh.vertices[triangle[0]][0])+','+str(self.Mesh.vertices[triangle[0]][1]))
                print(str(self.Mesh.vertices[triangle[1]][0])+','+str(self.Mesh.vertices[triangle[1]][1]))
                print(str(self.Mesh.vertices[triangle[2]][0])+','+str(self.Mesh.vertices[triangle[2]][1]))
            Vecx = ( (self.values[triangle[0]] - self.values[triangle[2]])*(self.Mesh.vertices[triangle[1]][1] - self.Mesh.vertices[triangle[0]][1])
                    + (self.values[triangle[1]] - self.values[triangle[0]])*(self.Mesh.vertices[triangle[2]][1] - self.Mesh.vertices[triangle[0]][1]) )/det
            Vecy = ( (self.values[triangle[0]] - self.values[triangle[1]])*(self.Mesh.vertices[triangle[2]][0] - self.Mesh.vertices[triangle[0]][0])
                    + (self.values[triangle[2]] - self.values[triangle[0]])*(self.Mesh.vertices[triangle[1]][0] - self.Mesh.vertices[triangle[0]][0]) )/det
            LTrianglesGrad.append([-Vecx/np.sqrt(Vecx**2 + Vecy**2),-Vecy/np.sqrt(Vecx**2 + Vecy**2)])
    return LTrianglesGrad

def showVectorField(self):
    L = self.computeGradientFlow()
    fig = go.Figure()
    self.Mesh.domain.addPlot(fig)
    fig.add_trace(go.Scatter(x=[P[0] for P in self.Mesh.vertices],
                            y = [P[1] for P in self.Mesh.vertices],
                    hoverinfo = "none",
                    showlegend = False,
                    mode="markers",
                    marker = dict(
                    color = self.values,
                    colorscale = 'viridis'
                    )))
    figQuiv = ff.create_quiver([(self.Mesh.vertices[T[0]][0]+self.Mesh.vertices[T[1]][0]+self.Mesh.vertices[T[2]][0])/3 for T in self.Mesh.triangles],
                                [(self.Mesh.vertices[T[0]][1]+self.Mesh.vertices[T[1]][1]+self.Mesh.vertices[T[2]][1])/3 for T in self.Mesh.triangles],
                                [V[0] for V in L], [V[1] for V in L])
    fig.add_traces(figQuiv.data)
    fig.update_layout(yaxis=dict(
        scaleanchor='x',
        scaleratio=1))
    fig.show()
