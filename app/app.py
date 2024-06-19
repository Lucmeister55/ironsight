from flask import Flask, request, render_template, redirect
import os

app = Flask(__name__)

@app.route('/', methods=['GET', 'POST'])
def index():
    if request.method == 'POST':
        return redirect('/progress')

    return render_template('index.html')

@app.route('/progress')  # Changed from bp.route to app.route
def progress():
    # Define the path to the progress.txt file
    templates_dir = os.path.join(app.root_path, 'templates')
    progress_file_path = os.path.join(templates_dir, 'progress.txt')

    # Check if the file exists
    if not os.path.exists(progress_file_path):
        return "No progress file found.", 404

    # Read the contents of the file
    with open(progress_file_path, 'r') as f:
        progress = f.read()

    # Return the contents of the file as a response
    return progress

@app.route('/results')
def results():
    return render_template('results.html')

if __name__ == '__main__':
    app.run(debug=True)