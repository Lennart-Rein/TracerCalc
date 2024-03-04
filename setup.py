from setuptools import setup, find_packages

VERSION = '0.0.1'
DESCRIPTION = 'Tracer planning'
LONG_DESCRIPTION = 'A lightweight library providing useful methods for planning tracer tests.'

# Setting up
setup(
    name="tracer_calc",
    version=VERSION,
    author="Lennart Rein",
    author_email="<lennartrein@gmail.com>",
    description=DESCRIPTION,
    long_description_content_type="text/markdown",
    long_description=LONG_DESCRIPTION,
    packages=find_packages(),
    install_requires=['matplotlib'],
    keywords=['python', 'tracer', 'groundwater', 'water', 'mass', "model"],
    classifiers=[
        "Development Status :: 1 - Planning",
        "Intended Audience :: Developers",
        "Programming Language :: Python :: 3",
        "Operating System :: Unix",
        "Operating System :: MacOS :: MacOS X",
        "Operating System :: Microsoft :: Windows",
    ]
)