from setuptools import setup, find_packages
import os

def read_version():
    version_file = os.path.join("bamurai", "VERSION")
    with open(version_file, encoding='utf-8') as f:
        return f.read().strip()

setup(
    name="bamurai",
    version=read_version(),
    packages=find_packages(),
    include_package_data=True,
    package_data={"bamurai": ["VERSION"]},
    entry_points={
        "console_scripts": [
            "bamurai=bamurai.cli:main",
        ],
    },
    license_files=('LICENSE',),
)
