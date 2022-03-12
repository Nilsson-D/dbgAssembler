from setuptools import setup, find_packages

setup(
    name="dbgAssembler",
    version="0.0.1",
    author="Daniel Nilsson",
    description="Simple genome assembler",
    packages=find_packages(),
    python_requires="3.9.10",
    install_requires=[
      "pygraphviz==1.9",
      "pandas==1.4.1",
      "numpy==1.22.2",
      "networkx==2.5",
      "matplotlib==3.5.1"
    ]

)
