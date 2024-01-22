from dash import Dash, dcc, html, Input, Output
import plotly.express as px
import pandas as pd
import numpy
import json

df = pd.read_csv("/home/franco/CLionProjects/untitled/sequential/points.csv", header=None)
neigh = pd.read_csv("/home/franco/CLionProjects/untitled/sequential/neighbours.csv", header=None)

""" fig = px.scatter_3d(x=df[0], y=df[1], z=df[2])
fig.show() """

app = Dash(__name__)

app.layout = html.Div([
    dcc.Graph(id="graph")
])


@app.callback(
    Output("graph", "figure"), 
    Input("graph", "clickData"))
def selectPoint(clickData):
    print(clickData)

    opacity = numpy.full(len(df), 0.5)

    if(clickData):
        pointId = clickData['points'][0]['pointNumber']
        neighbours = neigh.iloc[pointId].to_list()
        for n in neighbours:
            opacity[pointId] = 0.8
            opacity[n] = 1
        print(neighbours)
    fig = px.scatter_3d(x=df[0], y=df[1], z=df[2], color=opacity)

    return fig

app.run_server()
