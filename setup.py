"""
Setup script for selfies_analysis package
"""

from setuptools import setup, find_packages

with open("README.md", "r", encoding="utf-8") as fh:
    long_description = fh.read()

setup(
    name="selfies_analysis",
    version="0.1.0",
    author="Your Name",
    author_email="your.email@example.com",
    description="A package for analyzing SELFIES molecular representations",
    long_description=long_description,
    long_description_content_type="text/markdown",
    url="https://github.com/yourusername/selfies_analysis",
    packages=find_packages(),
    classifiers=[
        "Programming Language :: Python :: 3",
        "License :: OSI Approved :: MIT License",
        "Operating System :: OS Independent",
    ],
    python_requires=">=3.7",
    install_requires=[
        "numpy",
        "pandas",
        "matplotlib",
        "selfies",
        "rdkit",
        "scikit-learn",
    ],
)