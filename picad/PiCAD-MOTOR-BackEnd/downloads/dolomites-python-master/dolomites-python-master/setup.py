
from setuptools import setup


with open("README.md", "r", encoding="utf-8") as fh:
    long_description = fh.read()


setup(name='dolomites',
      version="0.0.1",
      author= "Luigi Alberti",
      author_email="luigi.alberti@unipd.it",
      description="A set of tools to design and simulate electric machine and drives",
      long_description=long_description,
      long_description_content_type="text/markdown",
      url="https://gitlab.com/LuigiAlberti/dolomites-python",
      packages=['dolomites'],
      classifiers=[
        "Programming Language :: Python :: 3",
        "License :: OSI Approved :: GNU General Public License v3 (GPLv3)",
        "Operating System :: OS Independent",
      ],
      install_requires=[
          'numpy',
          'pandas',
          'matplotlib',
          'scipy',
          'ezdxf',
          'PySide6'
          #'pytriangle @ git+https://git@github.com/pletzer/pytriangle.git@9c8537b2a2be94b7b85eb414930cc2c49a0f4368#egg=pytriangle'
          #'pytriangle',
          #'PySide2'
      ],
      python_requires='>=3.6'
      )
