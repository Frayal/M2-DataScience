# -*- coding: utf-8 -*-
import dash
import dash_core_components as dcc
import dash_html_components as html
import pandas as pd
import plotly.graph_objs as go
import numpy as np

app = dash.Dash()
df = pd.read_csv('scores.csv').iloc[1:]
if "istest" in df.columns.values:
    df['mois'] = df['istest'].apply(lambda x: '2017-12' if(x == 0) else  '2018-02' )
try:
    df1 = pd.read_csv('scores-test.csv').iloc[1:]
except Exception as e:
    df1 = pd.read_csv('scores.csv').iloc[1:]
try:
    df2 = pd.read_csv('scores-test1.csv').iloc[1:]
except Exception as e:
    df2 = pd.read_csv('scores.csv').iloc[1:]

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
    html.Label('Indicateur'),
        dcc.RadioItems(
        id='toggle1',
        options=[{'label': i, 'value': i} for i in ['Show', 'Hide']],
        value='Show'
        ),
        dcc.RadioItems(
        id='toggle2',
        options=[{'label': i, 'value': i} for i in ['Show', 'Hide']],
        value='Show'
        ),
        dcc.RadioItems(
        id='toggle3',
        options=[{'label': i, 'value': i} for i in ['Show', 'Hide']],
        value='Show'
        ),
        dcc.Dropdown(
        id= 'choice',
        options=[{'label':name,'value': name} for name in df.drop('mois',axis=1).columns.values]
            ,
        value = 'score TF1'
        ),
        dcc.Dropdown(
        id = 'date',
        options=[
            {'label': 'Decembre 2017', 'value': '2017-12'},
            {'label': 'Février 2018', 'value': '2018-02'},
            {'label': 'Mars 2018', 'value': '2018-03'}
        ],
        value=['2017-12','2018-02','2018-03'],
        multi=True
    ),
    html.Div(id='controls-container1', children=[
    dcc.Graph(
        id='score',
        figure={
            'data': [
                go.Scatter(
                    y=df[df['mois'] == '2017-12']['score TF1'],
                    mode='lines+markers',
                    name = '2017-12',
                    opacity=0.7,
                    marker={
                        'size': 15,
                        'line': {'width': 0.5, 'color': 'white'}
                    },
                ),
                go.Scatter(
                    y=df[df['mois'] == '2018-02']['score TF1'],
                    mode='lines+markers',
                    name = '2018-02',
                    opacity=0.7,
                    marker={
                        'size': 15,
                        'line': {'width': 0.5, 'color': 'white'}
                    },
                ),
                go.Scatter(
                    y=df[df['mois'] == '2018-03']['score TF1'],
                    mode='lines+markers',
                    name = '2018-03',
                    opacity=0.7,
                    marker={
                        'size': 15,
                        'line': {'width': 0.5, 'color': 'white'}
                    },
                )
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
    html.Div(id='controls-container2', children=[
    dcc.Graph(
        id='score1',
        figure={
            'data': [
                go.Scatter(
                    y=df1[df1['mois'] == '2017-12']['score TF1'],
                    mode='lines+markers',
                    name = 'train',
                    opacity=0.7,
                    marker={
                        'size': 15,
                        'line': {'width': 0.5, 'color': 'white'}
                    },
                ),
                go.Scatter(
                    y=df1[df1['mois'] == '2018-02']['score TF1'],
                    mode='lines+markers',
                    name = 'test',
                    opacity=0.7,
                    marker={
                        'size': 15,
                        'line': {'width': 0.5, 'color': 'white'}
                    },
                )
            ],
            'layout': go.Layout(
                xaxis={'title': 'Tours de l algorithme'},
                yaxis={'title': 'Score 1'},
                margin={'l': 40, 'b': 40, 't': 10, 'r': 10},
                legend={'x': 0, 'y': 1},
                hovermode='closest'
            )
        }
    ),]),
    html.Div(id='controls-container3', children=[
    dcc.Graph(
        id='score2',
        figure={
            'data': [
                go.Scatter(
                    y=df2[df2['mois'] == '2017-12']['score TF1'],
                    mode='lines+markers',
                    name = 'train',
                    opacity=0.7,
                    marker={
                        'size': 15,
                        'line': {'width': 0.5, 'color': 'white'}
                    },
                ),
                go.Scatter(
                    y=df2[df2['mois'] == '2018-02']['score TF1'],
                    mode='lines+markers',
                    name = 'test',
                    opacity=0.7,
                    marker={
                        'size': 15,
                        'line': {'width': 0.5, 'color': 'white'}
                    },
                )
            ],
            'layout': go.Layout(
                xaxis={'title': 'Tours de l algorithme'},
                yaxis={'title': 'Score 2'},
                margin={'l': 40, 'b': 40, 't': 10, 'r': 10},
                legend={'x': 0, 'y': 1},
                hovermode='closest'
            )
        }
    ) ]),
    ], style={'columnCount': 1})])

@app.callback(
    dash.dependencies.Output('score', 'figure'),
    [dash.dependencies.Input('choice','value'),
    dash.dependencies.Input('date','value')]
)

def Update(value,value_d):
    df = pd.read_csv('scores.csv').iloc[1:]
    print(len(value_d))
    return({'data': [
        go.Scatter(
            y=df[df['mois'] == str(v)][value],
            mode='lines+markers',
            name = str(v),
            opacity=0.7,
            marker={
                'size': 15,
                'line': {'width': 0.5, 'color': 'white'}
            },
        ) for v in value_d
    ],
    'layout': go.Layout(
        xaxis={'title': 'Tours de l algorithme'},
        yaxis={'title': 'Score'},
        margin={'l': 40, 'b': 40, 't': 10, 'r': 10},
        legend={'x': 0, 'y': 1},
        hovermode='closest'
    )
}
)
@app.callback(
    dash.dependencies.Output('score1', 'figure'),
    [dash.dependencies.Input('choice','value'),
    dash.dependencies.Input('date','value')]
)

def Update(value,value_d):
    print(len(value_d))
    return({'data': [
        go.Scatter(
            y=df1[df1['mois'] == str(v)][value],
            mode='lines+markers',
            name = str(v),
            opacity=0.7,
            marker={
                'size': 15,
                'line': {'width': 0.5, 'color': 'white'}
            },
        ) for v in value_d
    ],
    'layout': go.Layout(
        xaxis={'title': 'Tours de l algorithme'},
        yaxis={'title': 'Score 1'},
        margin={'l': 40, 'b': 40, 't': 10, 'r': 10},
        legend={'x': 0, 'y': 1},
        hovermode='closest'
    )
}
)
@app.callback(
    dash.dependencies.Output('score2', 'figure'),
    [dash.dependencies.Input('choice','value'),
    dash.dependencies.Input('date','value')]
)

def Update(value,value_d):
    print(len(value_d))
    return({'data': [
        go.Scatter(
            y=df2[df2['mois'] == str(v)][value],
            mode='lines+markers',
            name = str(v),
            opacity=0.7,
            marker={
                'size': 15,
                'line': {'width': 0.5, 'color': 'white'}
            },
        ) for v in value_d
    ],
    'layout': go.Layout(
        xaxis={'title': 'Tours de l algorithme'},
        yaxis={'title': 'Score 2'},
        margin={'l': 40, 'b': 40, 't': 10, 'r': 10},
        legend={'x': 0, 'y': 1},
        hovermode='closest'
    )
}
)
@app.callback(dash.dependencies.Output('controls-container1', 'style'), [dash.dependencies.Input('toggle1', 'value')])
def toggle_container(toggle_value):
    if toggle_value == 'Hide':
        return {'display': 'none'}
    else:
        return {'display': 'block'}
@app.callback(dash.dependencies.Output('controls-container2', 'style'), [dash.dependencies.Input('toggle2', 'value')])
def toggle_container(toggle_value):
    if toggle_value == 'Hide':
        return {'display': 'none'}
    else:
        return {'display': 'block'}
@app.callback(dash.dependencies.Output('controls-container3', 'style'), [dash.dependencies.Input('toggle3', 'value')])
def toggle_container(toggle_value):
    if toggle_value == 'Hide':
        return {'display': 'none'}
    else:
        return {'display': 'block'}

if __name__ == '__main__':
    app.run_server(debug=True,port=3003)
