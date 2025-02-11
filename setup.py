from setuptools import setup, find_packages

def read_version():
    with open("VERSION") as version_file:
        return version_file.read().strip()

setup(
    name="bamurai",
    version=read_version(),
    packages=find_packages(),
    include_package_data=False,
    entry_points={
        "console_scripts": [
            "bamurai=bamurai.cli:main",
        ],
    },
    license_files=('LICENSE',),
)