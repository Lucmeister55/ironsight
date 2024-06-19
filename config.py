import os
import pkg_resources

class Config:
    SECRET_KEY = os.urandom(24)
    PKG_FOLDER = pkg_resources.get_distribution("ironsight").location
    TEMP_FOLDER = os.path.join(PKG_FOLDER, 'app/temp')
    STATIC_FOLDER = os.path.join(PKG_FOLDER, 'app/static')
    RESULTS_FOLDER = os.path.join(PKG_FOLDER, 'app/static/results')
    TEMPLATES_FOLDER = os.path.join(PKG_FOLDER, 'app/templates')

config = Config()
