import dash
import dash_cytoscape as cyto
import dash_html_components as html

app = dash.Dash(__name__)
app.layout = html.Div([
    cyto.Cytoscape(
        id='cytoscape',
        elements=[
            {'data': {'id': 'one', 'label': 'Node 1'}, 'position': {'x': 50, 'y': 50}},
            {'data': {'id': 'two', 'label': 'Node 2'}, 'position': {'x': 200, 'y': 200}},
            {'data': {'source': 'one', 'target': 'two','label': 'Node 1 to 2'}},
            {'data': {'source': 'one', 'target': 'one','label': 'Node 1 to 1'}},
            {'data': {'source': 'two', 'target': 'one','label': 'Node 2 to 1'}}
        ],
        layout={'name': 'preset'}
    )
    cy.edges('edge').style({
    "curve-style": "bezier",
    "target-arrow-shape": "triangle"
})
])

if __name__ == '__main__':
    app.run_server(debug=False)


# import dash
#
# from demos.editor.callbacks import assign_callbacks
# from demos.editor.layout import layout as cytoscape_layout
#
# app = dash.Dash(__name__)
# server = app.server
#
# app.scripts.config.serve_locally = True
# app.css.config.serve_locally = True
#
#
# app.layout = cytoscape_layout
# assign_callbacks(app)
#
#
# if __name__ == '__main__':
#     app.run_server(debug=True)
