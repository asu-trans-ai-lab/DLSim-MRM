import setuptools

with open("README.md", "r", encoding="utf-8") as fh:
     long_description = fh.read()

setuptools.setup(
    name="DLSim", 
    version="0.2.1",
    author="Dr.Cafer Avci, Zhiqiang Niu, Jiawei Liu, Dr.Xuesong Zhou",
    author_email="author@example.com",
    License="Apache Software License (http://www.apache.org/licenses/LICENSE-2.0)",
    description="demo",
    long_description=long_description,
    long_description_content_type="text/markdown",
    url="https://github.com/asu-trans-ai-lab/DLSim",
    packages=['DLSim'],
    install_requires=[            
          'numpy',
    ],
    classifiers=[
        "License :: OSI Approved :: Apache Software License",
        "Operating System :: OS Independent",
        "Programming Language :: Python :: 3",
        "Programming Language :: Python :: Implementation :: PyPy",
    ],
)