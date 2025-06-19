"""
Setup script for selfies_analysis package
"""

from setuptools import setup, find_packages
import os

# Read the contents of README file
this_directory = os.path.abspath(os.path.dirname(__file__))
with open(os.path.join(this_directory, "README.md"), "r", encoding="utf-8") as fh:
    long_description = fh.read()

# Base requirements that can be installed via pip
base_requirements = [
    "numpy>=1.19.0",
    "pandas>=1.1.0", 
    "matplotlib>=3.3.0",
    "selfies>=2.0.0",
    "scikit-learn>=0.24.0",
]

# RDKit can be tricky - make it optional for now
rdkit_requirements = [
    "rdkit-pypi>=2022.9.1",  # This is the pip-installable version
]

setup(
    name="selfies-analysis",
    version="0.4.1",
    author="Your Name",  # TODO: Replace with actual name
    author_email="your.email@example.com",  # TODO: Replace with actual email
    description="A package for analyzing SELFIES molecular representations",
    long_description=long_description,
    long_description_content_type="text/markdown",
    url="https://github.com/yourusername/selfies_analysis",  # TODO: Replace with actual URL
    packages=find_packages(),
    classifiers=[
        "Development Status :: 4 - Beta",
        "Intended Audience :: Science/Research",
        "Topic :: Scientific/Engineering :: Chemistry",
        "Topic :: Scientific/Engineering :: Artificial Intelligence",
        "Programming Language :: Python :: 3",
        "Programming Language :: Python :: 3.7",
        "Programming Language :: Python :: 3.8", 
        "Programming Language :: Python :: 3.9",
        "Programming Language :: Python :: 3.10",
        "Programming Language :: Python :: 3.11",
        "License :: OSI Approved :: MIT License",
        "Operating System :: OS Independent",
    ],
    python_requires=">=3.7",
    install_requires=base_requirements,
    extras_require={
        "rdkit": rdkit_requirements,
        "all": base_requirements + rdkit_requirements,
        "dev": [
            "pytest>=6.0",
            "pytest-cov>=2.0",
            "black>=21.0",
            "flake8>=3.8",
            "mypy>=0.800",
        ],
    },
    keywords="chemistry, molecules, selfies, smiles, rdkit, cheminformatics, analysis",
    project_urls={
        "Bug Reports": "https://github.com/yourusername/selfies_analysis/issues",
        "Source": "https://github.com/yourusername/selfies_analysis",
        "Documentation": "https://github.com/yourusername/selfies_analysis#readme",
    },
    include_package_data=True,
)