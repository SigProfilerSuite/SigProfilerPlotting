from setuptools import setup

VERSION_TEMPLATE = "version = '{version}'\n"


with open("README.md") as f:
    readme = f.read()

setup(
    name="sigProfilerPlotting",
    use_scm_version={
        "write_to": "sigProfilerPlotting/_version.py",
        "write_to_template": VERSION_TEMPLATE,
        "fallback_version": "0+unknown",
    },
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
        "sigProfilerPlotting.sbs4608",
    ],
    python_requires=">=3.9",
    install_requires=[
        "matplotlib>=3.4.3",
        "pandas>=2.0.0",
        "scikit-learn>=1.1.3",
        "pillow>=10.0.0",
        "plotly>=6.1.1",
        "kaleido>=1.0.0",
        "pypdf>=5.0.0",
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
