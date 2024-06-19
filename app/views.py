from flask import current_app as app
from flask import Blueprint, render_template, request, redirect, url_for
from ironsight import Predict_Ferro_Betas
from threading import Thread
import os
from flask import session
from concurrent.futures import ThreadPoolExecutor
from uuid import uuid4
import pandas as pd
import plotly.graph_objects as go

executor = ThreadPoolExecutor(max_workers=1)
futures = {}

bp = Blueprint('main', __name__)

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


def run_prediction(file_path, id, model_type, progress_lines_path, progress_percentage_path, results_dir, temp_dir, pkg_dir, output_html):

    # Call the processing function
    results_csv_path, explainer = Predict_Ferro_Betas.predict([file_path], [id], model_type, progress_lines_path, progress_percentage_path, results_dir, temp_dir, pkg_dir)

    create_results_html(results_csv_path, output_html, explainer)
    # Set a flag in the session to indicate that the prediction is done
    return True

@bp.route('/', methods=['GET', 'POST'])
def index():
    if request.method == 'POST':
        file_path = request.form['file_path']
        id = request.form['id']
        model_type = request.form['model_type']
        outdir = request.form['outdir']

        # Check if the file exists
        if not os.path.exists(file_path):
            return f"File {file_path} does not exist.", 400

        results_dir = app.config['RESULTS_FOLDER']
        temp_dir = app.config['TEMP_FOLDER']
        pkg_dir = app.config['PKG_FOLDER']

        templates_dir = os.path.join(app.root_path, 'templates')
        progress_lines_path = os.path.join(templates_dir, 'progress_lines.txt')
        progress_percentage_path = os.path.join(templates_dir, 'progress_percentage.txt')
        output_html = os.path.join(app.config['TEMPLATES_FOLDER'], 'results.html')

        with open(progress_lines_path, 'w') as progress_file:
            progress_file.write('Starting...' + '\n')

        with open(progress_percentage_path, 'w') as progress_file:
            progress_file.write('0' + '\n')

        # Start the prediction process in a new thread and get the Future object
        future = executor.submit(run_prediction, file_path, id, model_type, progress_lines_path, progress_percentage_path, results_dir, temp_dir, pkg_dir, output_html)

        # Generate a unique task ID
        task_id = str(uuid4())

        # Store the Future object in the global dictionary
        futures[task_id] = future

        # Store the task ID in the session
        session['task_id'] = task_id

        return redirect(url_for('main.progress'))

    return render_template('index.html')

@bp.route('/progress')
def progress():
    templates_dir = os.path.join(app.root_path, 'templates')
    progress_lines_path = os.path.join(templates_dir, 'progress_lines.txt')
    progress_percentage_path = os.path.join(templates_dir, 'progress_percentage.txt')

    with open(progress_lines_path, 'r') as f:
        progress_lines = f.read()

    with open(progress_percentage_path, 'r') as f:
        progress_percentage = f.read()

    # Replace newlines with <br> for HTML
    progress_lines = progress_lines.replace('\n', '<br>')

    # Get the task ID from the session
    task_id = session.get('task_id')

    # Get the Future object from the global dictionary
    future = futures.get(task_id)

    # Check if the future is done
    prediction_done = future.done() if future else False

    return render_template('progress.html', progress_lines=progress_lines, progress_percentage=progress_percentage, prediction_done=prediction_done)

@bp.route('/results')
def results():
    return render_template('results.html')