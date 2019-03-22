# -*- coding: utf-8 -*-
import dash
import dash_core_components as dcc
import dash_html_components as html
import pandas as pd
import plotly.graph_objs as go

def create_app():
    app = dash.Dash(__name__)

    df = pd.read_csv('gdp-life-exp-2007.csv')
    app.layout = html.Div(children=[
        html.H1(children='Dash works on DNAnexus!'),

        dcc.Graph(
            id='life-exp-vs-gdp',
            figure={
                'data': [
                    go.Scatter(
                        x=df[df['continent'] == i]['gdp per capita'],
                        y=df[df['continent'] == i]['life expectancy'],
                        text=df[df['continent'] == i]['country'],
                        mode='markers',
                        opacity=0.7,
                        marker={
                            'size': 15,
                            'line': {'width': 0.5, 'color': 'white'}
                        },
                        name=i
                    ) for i in df.continent.unique()
                ],
                'layout': go.Layout(
                    xaxis={'type': 'log', 'title': 'GDP Per Capita'},
                    yaxis={'title': 'Life Expectancy'},
                    legend={'x': 0, 'y': 1},
                    hovermode='closest',
                    title="As proof, here is an interactive Gapminder-style scatter plot"
                )
            }
        ),
        dcc.Markdown('''
You can also write in Markdown, so we can easily write documentation straight into the interface. This is how you make an applet open up HTTPS by the way. Just add this to the dxapp.json:
```
"httpsApp": {"ports":[443], "shared_access": "VIEW"},
```
And then your web app should output on port 443.
        '''),
        dcc.Markdown('''
For more information on what you can build with Dash, see the Dash [tutorial](https://dash.plot.ly/).
        ''')
    ])
    return app
