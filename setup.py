import os
import shutil

from setuptools import setup
from setuptools.command.install import install

# remove the dist folder first if exists
if os.path.exists("dist"):
    shutil.rmtree("dist")


def readme():
    with open("README.rst") as f:
        return f.read()


VERSION = "1.4.3"


def write_version_py(filename="sigProfilerPlotting/version.py"):
    # Copied from numpy setup.py
    cnt = """
# THIS FILE IS GENERATED FROM SIGPROFILERPLOTTING SETUP.PY
short_version = '%(version)s'
version = '%(version)s'
update = 'v1.4.3: Fix pandas 3.12 compatibility issues in tmbplot.py'
    
    """
    fh = open(filename, "w")
    fh.write(
        cnt
        % {
            "version": VERSION,
        }
    )
    fh.close()


write_version_py()


with open("README.md") as f:
    readme = f.read()

setup(
    name="sigProfilerPlotting",
    version=VERSION,
    description="SigProfiler plotting tool",
    long_description=readme,
    long_description_content_type="text/markdown",
    url="https://github.com/alexandrovlab/SigProfilerPlotting",
    author="Erik Bergstrom",
    author_email="ebergstr@eng.ucsd.edu",
    license="UCSD",
    packages=[
        "sigProfilerPlotting",
        "sigProfilerPlotting.reference_formats",
        "sigProfilerPlotting.fonts",
        "sigProfilerPlotting.controllers",
    ],
    python_requires=">=3.9",
    install_requires=[
        "matplotlib>=3.4.3",
        "pandas>=2.0.0",
        "scikit-learn>=1.1.3",
        "pillow>=10.0.0",
    ],
    extras_require={
        "tests": [
            "pytest",
            "scikit-image>=0.21.0",
            "numpy>=2.0.0",
        ],
    },
    entry_points={
        "console_scripts": [
            "SigProfilerPlotting=sigProfilerPlotting.sigProfilerPlotting_CLI:main_function",
        ],
    },
    package_data={"": ["fonts/*.ttf"]},
    include_package_data=True,
    zip_safe=False,
)
