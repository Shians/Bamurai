from setuptools import setup, find_packages

setup(
    name="bamurai",
    version="0.1.0",
    packages=find_packages(),
    include_package_data=False,
    entry_points={
        "console_scripts": [
            "bamurai=bamurai.cli:main",
        ],
    },
    license_files=('LICENSE',),
)