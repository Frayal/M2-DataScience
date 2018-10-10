# -*- coding: utf-8 -*-
import dash
import dash_core_components as dcc
import dash_html_components as html
import pandas as pd
import plotly.graph_objs as go
import numpy as np
from datetime import datetime as dt
from datetime import date, timedelta

app = dash.Dash()
df = pd.read_csv('cout.csv')


colors = {
    'background': '#0000',
    'text': '#7FDBFF'
}

app.layout = html.Div(style={'backgroundColor': colors['background']}, children=[
    html.H1(
        children='Bonjour Alexis',
        style={
            'textAlign': 'center',
            'color': colors['text']
        }
    ),

    html.Div(children='Visualisation des résultats', style={
        'textAlign': 'center',
        'color': colors['text']
    }),
    html.Div([
    html.Label('Date'),
        dcc.DatePickerRange(
        id='date-input',
        min_date_allowed=dt(2017, 12, 1),
        max_date_allowed=dt(2018, 4, 1),
        initial_visible_month=dt(2017, 12, 1),
        end_date=dt(2017,12,31)
        ),
    html.Label('Chaîne'),
    dcc.Dropdown(
        id = 'Chaîne',
        options = [{'label': i, 'value': i} for i in ['TF1','M6','France 2','France 3']],
        value=['TF1'],
        multi=True
    ),
    html.Label('Part of the day'),
        dcc.RadioItems(
        id='toggle1',
        options=[{'label': i, 'value': i} for i in ['tout','matinee','apresmidi','soiree']],
        value='tout'
    ),
    html.Div(id='controls-container1', children=[
    dcc.Graph(
        id='score',
        figure={
            'data': [
                go.Scatter(
                    y = df['2017-12-02_TF1_tout'],
                    mode='lines+markers',
                    name = 'TF1',
                    opacity=0.7,
                    marker={
                        'size': 15,
                        'line': {'width': 0.5, 'color': 'white'}
                    },
                ),

            ],
            'layout': go.Layout(
                xaxis={'title': 'Tours de l algorithme'},
                yaxis={'title': 'Score'},
                margin={'l': 40, 'b': 40, 't': 10, 'r': 10},
                legend={'x': 0, 'y': 1},
                hovermode='closest'
            )
        }
    ),]),

    ], style={'columnCount': 1})])


@app.callback(
    dash.dependencies.Output('controls-container1', 'children'),
    [dash.dependencies.Input('date-input','start_date'),
    dash.dependencies.Input('date-input','end_date'),
    dash.dependencies.Input('Chaîne','value'),
    dash.dependencies.Input('toggle1','value')]
)

def Update(start,end,value_d,value_c):
    df = pd.read_csv('cout.csv')
    a = date(int(start.split('-')[0]),int(start.split('-')[1]),int(start.split('-')[2]))
    b = date(int(end.split('-')[0]),int(end.split('-')[1]),int(end.split('-')[2]))
    l = []
    for i in range((b-a).days+1):
        here = []
        not_here = []
        for v in value_d:
            try:
                y = df[str(a+timedelta(i))+'_'+v+'_'+value_c]
                here.append(v)
            except Exception as e:
                not_here.append(v)
        if(len(here)==0):
            pass
        else:
            l.append(dcc.Graph(
                id='score'+str(i),
                figure={
                    'data': [
                    go.Scatter(
                        y=df[str(a+timedelta(i))+'_'+v+'_'+value_c],
                        mode='lines+markers',
                        name = str(v),
                        opacity=0.7,
                        marker={
                            'size': 15,
                            'line': {'width': 0.5, 'color': 'white'}
                        },
                    ) for v in here],
                    'layout': go.Layout(
                        xaxis={'title': 'Tours de l algorithme pour le '+str(a+timedelta(i))},
                        yaxis={'title': 'Score'},
                        margin={'l': 40, 'b': 40, 't': 10, 'r': 10},
                        legend={'x': 0, 'y': 1},
                        hovermode='closest'
                    )
                }
            ))

    return(l)

if __name__ == '__main__':
    app.run_server(debug=True,port=3004)
