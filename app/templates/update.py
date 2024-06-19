import pandas as pd
import plotly.graph_objects as go

def create_results_html(file_path, output_path, explainer):
    # Load the data from CSV
    df = pd.read_csv(file_path)

    # Create traces
    trace1 = go.Bar(y=df['id'], x=df['probability_S'], name='Sensitive', orientation='h', marker=dict(color='green'), hovertemplate = '%{x:.2f}')
    trace2 = go.Bar(y=df['id'], x=df['probability_R'], name='Resistant', orientation='h', marker=dict(color='red'), hovertemplate = '%{x:.2f}')

    # Create the layout
    layout = go.Layout(
        barmode='stack', 
        title='MODEL PREDICTION',
        xaxis_title='Probability',
        yaxis_title='Sample ID',
        plot_bgcolor='rgba(0,0,0,0)',  # make the plot background transparent
        xaxis=dict(showgrid=False),  # remove the x-axis grid lines
        yaxis=dict(showgrid=False),  # remove the y-axis grid lines
        height=300,  # adjust the height of the plot
        width=600,  # adjust the width of the plot
    )

    # Create the figure and add traces
    fig = go.Figure(data=[trace1, trace2], layout=layout)

    # Save the plot as an HTML file
    fig.write_html(output_path)

file_path = "/data/lvisser/ironsight/app/static/results/results.csv"
output_path = "/data/lvisser/ironsight/app/templates/results.html"

create_results_html(file_path, output_path, explainer)