import setuptools


setuptools.setup(
    name="purplebook",
    version="0.1.0",
    author="gully",
    author_email="igully@gmail.com",
    description="Scripts and utilities from Princples of Planetary Climate",
    long_description="A Python 3 port of Scripts and utilities from Princples of Planetary Climate",
    long_description_content_type="text/markdown",
    url="https://github.com/BrownDwarf/purplebook",
    install_requires=[
        "numpy",
        "scipy",
        "astropy>=4.1",
        "specutils",
        "pandas",
        "matplotlib",
    ],
    packages=setuptools.find_packages(where="src"),
    package_dir={"": "src"},
    package_data={
        # If any package contains *.csv files, include them:
        "": ["*.csv"]
    },
    classifiers=[
        "Programming Language :: Python :: 3",
        "License :: OSI Approved :: MIT License",
        "Operating System :: OS Independent",
    ],
    python_requires=">=3.7",
)
