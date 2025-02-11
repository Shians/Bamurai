import os

def get_version():
    version_file = os.path.join(os.path.dirname(__file__), "VERSION")
    with open(version_file) as f:
        return f.read().strip()
