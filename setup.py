from setuptools import setup, find_packages

setup(
    name="autopocket2crest",
    version="1.0.0",
    author="Christian Fellinger",
    author_email="christian.fellinger@univie.ac.at",
    description="Automated protein–ligand pocket preparation pipeline for CREST conformational sampling.",
    packages=find_packages(),
    entry_points={
        "console_scripts": ["autopocket2crest=autopocket2crest.cli:main"]
    },
    python_requires=">=3.8",
)
